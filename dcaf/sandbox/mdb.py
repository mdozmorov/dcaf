import collections
import struct

class Struct(type):
    @classmethod
    def __prepare__(metacls, name, bases):
        return collections.OrderedDict()

    def __new__(cls, name, bases, classdict):
        c = type.__new__(cls, name, bases, classdict)
        #c._orderedKeys = classdict.keys()

        format = "="
        keys = []
        for k,v in classdict.items():
            if k.startswith("_"):
                continue
            keys.append(k)
            format += v

        unpacker = struct.Struct(format)
        
        def __init__(self, buffer, offset=0):
            self.size = unpacker.size
            data = unpacker.unpack_from(buffer, offset)
            for k,v in zip(keys, data):
                setattr(self, k, v)

        c.__init__ = __init__
        return c

from datetime import datetime, timedelta
import mmap
import struct
import enum
import io
import array
import uuid

from collections import namedtuple

MAX_OBJ_NAME_SIZE = 256

class PageType(enum.Enum):
    DBDefinition = 0x00
    Data = 0x01
    TableDefinition = 0x02
    IntermediateIndex = 0x03
    LeafIndex = 0x04
    PageUsage = 0x05

OLE_TIME_ZERO = datetime(1899, 12, 30, 0, 0, 0)

class ColumnType(enum.Enum):
    BOOL = (0x01, 1)
    BYTE = (0x02, 8)
    INT = (0x03, 16, lambda x: int.from_bytes(x, "little"))
    LONGINT = (0x04, 32, lambda x: int.from_bytes(x, "little"))
    MONEY = (0x05, 64)
    FLOAT = (0x06, 32, lambda x: struct.unpack("<f", x)[0])
    DOUBLE = (0x07, 64, lambda x: struct.unpack("<d", x)[0])
    # FIXME: Figure out how to convert DATETIME
    # Supposedly, it is a double which represents
    # time in days since OLE_TIME_ZERO
    DATETIME = (0x08, 64, lambda x: OLE_TIME_ZERO)# + \
                #timedelta(days=struct.unpack("<d", x)[0]))
    BINARY = (0x09, 255 * 8) # or 255?
    TEXT = (0x0A, 255 * 8)
            #lambda x: x.decode("utf16"))
            #lambda x: x[:x.find(b"\x00")])
    OLE = (0x0B, 0)
    MEMO = (0x0C, 0)
    UNKNOWN_0D = (0x0D, 0)
    UNKNOWN_0E = (0x0E, 0)
    GUID = (0x0F, 128, lambda x: uuid.UUID(bytes=x))
    NUMERIC = (0x10, 17)

    def __init__(self, flag, size, ctor=lambda x: x):
        self.flag = flag
        self.size = size
        self.ctor = ctor
        
    @property
    def fixed(self):
        return self.size > 0
    
class TableDef(metaclass=Struct):
    tdef_len = "I"
    unknown1 = "I"
    num_rows = "I"
    autonumber = "I"
    autonum_flag = "B"
    unknown2 = "3s"
    ct_autonum = "I"
    unknown3 = "Q"
    table_type = "B"
    max_cols = "H"
    num_var_cols = "H"
    num_cols = "H"
    num_idx = "I"
    num_real_idx = "I"
    used_pages = "I"
    free_pages = "I"

class ColumnDef(metaclass=Struct):
    col_type = "B"
    unknown1 = "I"
    col_num = "H" # Includes deleted columns
    offset_V = "H"
    col_num2 = "H" # Doesn't include deleted columns
    misc = "H"
    misc_ext = "H"
    bitmask = "b"
    misc_flags = "b"
    unknown2 = "I"
    offset_F = "H"
    col_len = "H"
    
class DataPageDef(metaclass=Struct):
    page_type = "B"
    unknown1 = "B"
    free_space = "H"
    tdef_pg = "I"
    unknown2 = "I" # Only for JET4
    num_rows = "H"

MDBConstants = namedtuple("MDBConstants", """
pg_size
row_count_offset 
tab_num_rows_offset
tab_num_cols_offset
tab_num_idxs_offset
tab_num_ridxs_offset
tab_usage_map_offset
tab_first_dpg_offset
tab_cols_start_offset
tab_ridx_entry_size
col_fixed_offset
col_size_offset
col_num_offset
tab_col_entry_size
tab_free_map_offset
tab_col_offset_var
tab_col_offset_fixed
tab_row_col_num_offset""")

MDB4Constants = MDBConstants(4096, 0x0c, 16, 45, 
                             47, 51, 55, 56, 63, 
                             12, 15, 23, 5, 25, 
                             59, 7, 21, 9)
MDB3Constants = MDBConstants(2048, 0x08, 12, 25, 
                             27, 31, 35, 36, 43, 
                             8, 13, 16, 1, 18, 39, 
                             3, 14, 5)

class MDB(object):
    def __init__(self, path):
        self._handle = io.open(path, "r+b")
        self._map = mmap.mmap(self._handle.fileno(), 0)
        
        self._map.seek(0x14)
        self.version = int.from_bytes(self._map.read(4), 
                                      "little") + 3
        self._constants = MDB4Constants \
                          if self.version == 4 \
                             else MDB3Constants
        self._endianness = "little"
        self._catalog = Table(self, 2)
    
    def _read(self, offset, format="=I"):
        self._map.seek(offset)
        st = struct.Struct(format)
        return st.unpack(self._map.read(st.size))[0]

    def table(self, name):
        pass

class Table(object):
    def __init__(self, mdb, page):
        self._mdb = mdb

        fmt = mdb._constants
        start = pos = (fmt.pg_size * page)
        # Required for valid table def page
        assert mdb._map[pos] == 0x02

        pos += 8
        self._def = TableDef(mdb._map, offset=pos)

        pos += self._def.size + \
              (self._def.num_real_idx * \
               fmt.tab_ridx_entry_size)

        ncol = self._def.num_cols

        self._columns = []
        for i in range(ncol):
            col = ColumnDef(mdb._map, offset=pos)
            for col_type in ColumnType:
                if col_type.flag == col.col_type:
                    col.col_type = col_type
            self._columns.append(col)
            pos += fmt.tab_col_entry_size

        mdb._map.seek(pos)
        for i,col in enumerate(self._columns):
            if mdb.version == 4:
                name_sz = struct.unpack("=H", 
                                        mdb._map.read(2))[0]
            else:
                raise NotImplementedError
            assert name_sz <= MAX_OBJ_NAME_SIZE
            name = mdb._map.read(name_sz)
            col.name = name.decode("utf16") \
                   if mdb.version == 4 \
                      else name.decode("ascii")

        self._columns.sort(key=lambda x: x.col_num)
        self._data_pg = struct.unpack_from("=H", mdb._map,
                                           start+fmt.tab_first_dpg_offset)[0]

    
    def __iter__(self):
        page = self._data_pg
        # FIXME: iterate over multiple pages
        for col in self._columns:
            print(col.name, col.col_type, col.bitmask & 0x01, col.offset_V, col.offset_F, col.col_len)
        return self._iter_page(page)
        #return self._iter_page(14)
   
    def _iter_page(self, page):
        mm = self._mdb._map
        pos = base = page * self._mdb._constants.pg_size
        pdef = DataPageDef(mm, pos)
        pos += pdef.size
        
        data = mm[pos:(pos+pdef.num_rows*2)]
        offsets = list(reversed(array.array("H", data)))
        #offsets.append(self._mdb._constants.pg_size - 1)
        offsets.append(self._mdb._constants.pg_size)
        print(offsets)
        
        for start, end in zip(offsets[:-1], offsets[1:]):
            deleted = 0x8000 & start
            pointer = 0x4000 & start
            if pointer or deleted:
                # FIXME: handle pointers
                continue
            start = (start & 0x1fff) + base
            end = (end & 0x1fff) + base

            row = [None] * len(self._columns)
            num_cols = struct.unpack_from("=H", mm, 
                                          start)[0]

            # FIXME: a row has different
            # structure when it has no variable columns

            null_mask_len = (num_cols + 7) // 8
            pos = end - null_mask_len - 2
            var_len = struct.unpack_from("=H", mm, pos)[0]
            null_mask = mm[pos:(pos+null_mask_len)]
            print(start, end, num_cols, 
                  null_mask_len, var_len)

            # Read the variable columns
            if var_len:
                var_offset_data = mm[(pos-(var_len*2)-2):pos]
                var_offsets = struct.unpack("="+"H"*(var_len+1), 
                                            var_offset_data)
                var_offsets = [start+o for o in
                               reversed(var_offsets)]
                var_ix = 0
                for col in self._columns:
                    type = col.col_type
                    if not type.fixed:
                        v_start = var_offsets[var_ix]
                        v_end = var_offsets[var_ix+1]
                        if v_start < v_end:
                            data = mm[v_start:v_end]
                            if type == ColumnType.MEMO and len(data) >= 12:
                                memo_len, lval_dp, unknown = struct.unpack("=III", data[:12])
                                bitmask = memo_len >> 24
                                memo_len &= 0x00ffffff
                                if 0x80 & bitmask:
                                # Otherwise it requires an LVAL page
                                # which is a pain
                                    row[col.col_num] = data[12:]
                            else:
                                row[col.col_num] = type.ctor(data)
                        # Otherwise, it'll remain None
                        var_ix += 1
            else:
                pass

            for col in self._columns:
                type = col.col_type
                if col.bitmask & 0x01:
                    mm.seek(start+col.offset_F)
                    data = mm.read(col.col_len)
                    row[col.col_num] = type.ctor(data)
                if type == ColumnType.BOOL:
                    pass
            yield row

path = "/home/gilesc/data/SORD29.mdb"
db = MDB(path)
for i,row in enumerate(db._catalog):
    print(row[:3])
    if i == 5:
        break
