#pragma once

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <string>
#include <iostream>

namespace dcaf {

  class MMFile {
    int handle;
    off_t fileSize;
    void *data;

  public:
    MMFile(std::string path) {
      struct stat st;
      stat(path.c_str(), &st);
      this->fileSize = st.st_size;

      handle = open(path.c_str(), O_RDONLY);
      this->data = mmap(0, fileSize, PROT_READ, 
                        MAP_SHARED, handle, 0);
      if (this->data == MAP_FAILED) {
        throw 5;
      }
    }

    ~MMFile() {
      munmap(this->data, this->fileSize);
      close(handle);
    }

    void* operator[](size_t offset) {
      return ((char*) data) + offset;
    }
  };
}
