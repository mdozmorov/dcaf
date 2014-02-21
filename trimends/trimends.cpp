#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>

using namespace std;

float MISMATCH_RATIO = 0.05;
int MIN_DUPLICATES = 5;
int MIN_LENGTH = 20;
bool VERBOSE = false;

struct FASTQ {
    string id;
    string seq;
    string qual;
};

FASTQ* read_fastq(istream& strm) {
    string ignore;
    FASTQ* result = new FASTQ;
    if (!getline(strm, result->id))
        return NULL;
    int lpos = result->id.find(" length"); 
    if (lpos == string::npos)
        lpos = result->id.length();
    result->id = result->id.substr(1,lpos);
    getline(strm, result->seq);
    getline(strm, ignore);
    getline(strm, result->qual);
    return result;
}

int trim(FASTQ* fq) {
    int i=1;
    char b = fq->seq[0];
    float nb = 1;
    float nnb = 0;
    float ratio = 0;
    while (ratio < MISMATCH_RATIO) {
        if (fq->seq[i]==b)
            nb++;
        else
            nnb++;
        ratio = nnb / nb;
        ++i;
    }
    return i > MIN_DUPLICATES ? i-1 : 0;
}

void usage() {
    cerr << "\nUSAGE:\n./trimends [options] fastq-file\n\n";
    cerr << "-v Verbose. Print any reads that were modified or discarded to stderr.\n";
    cerr << "-n [int] Min duplicated bases on an end to be considered for clipping.\n";
    cerr << "-r [0-100] Max ratio of mismatched bases to be considered for clipping.\n";
    cerr << "-l [int] Minimum length of a post-clipped read to avoid discarding.\n";
    cerr << "\n\n";
    exit(1);
}

int main(int argc, char* argv[]) {
    if (argc == 1)
        usage();

    char c;
    while ((c = getopt(argc,argv,"vr:n:l:")) != -1) {
        switch (c) {
            case 'v':
                VERBOSE = true;
                break;
            case 'r':
                MISMATCH_RATIO = atof(optarg) / 100.0;
                if (MISMATCH_RATIO < 0 || MISMATCH_RATIO > 1) {
                    cerr << "Invalid argument: mismatch ratio (-r) must be between 0 and 100\n";
                    exit(1);
                }
                break;
            case 'n':
                MIN_DUPLICATES = atoi(optarg);
                break;
            case 'l':
                MIN_LENGTH = atoi(optarg);
                break;
            case 'h':
                usage();
            default:
                usage();
        }
    }

    //if (optind >= argc)
     //   usage();

    //ifstream strm(argv[optind]);
    FASTQ* fq;
    map<int, int> clipped_counts;

    int nreads = 0;
    int nclipped1side = 0;
    int nclipped2side = 0;
    int nrepetitive = 0;
    int nbelowmin = 0;

    while (fq = read_fastq(cin)) {
        nreads++;
        string seq = fq->seq;
        int t = trim(fq);
        reverse(fq->seq.begin(), fq->seq.end());
        int rt = fq->seq.length() - trim(fq);
        if (rt < t) {
            if (VERBOSE) 
                cerr << "REPETITIVE: \n" << fq->seq << "\n\n";
            nrepetitive++;
            continue;
        }
        reverse(fq->seq.begin(), fq->seq.end());
        int len;
        if (t || rt-fq->seq.length()) {
            int nclipped = fq->seq.length() - (rt - t);
            clipped_counts[nclipped]++;
            if (VERBOSE) {
                cerr << fq->seq << endl;
            }
            fq->seq = fq->seq.substr(t, rt - t);
            fq->qual = fq->qual.substr(t, rt - t);
            len = fq->seq.length();
            if (VERBOSE) {
                cerr << fq->seq << endl << endl;
            }
            if (len >= MIN_LENGTH) {
                if (t && rt)
                    nclipped2side++;
                else
                    nclipped1side++;
            }
        }  else {
            len = fq->seq.length();
        }
        if (len < MIN_LENGTH) {
            if (VERBOSE)
                cerr << "TOO SHORT:\n" << fq->seq << endl << endl;
            nbelowmin++;
            continue;
        }
        cout << "@" << fq->id << endl; //" length=" << len << endl;
        cout << fq->seq << endl;
        cout << "+\n"; // << fq->id << " length=" << len << endl;
        cout << fq->qual << endl; 
        delete fq;
    }

    cerr << nreads << " reads were processed.\n";
    cerr << nrepetitive << " reads were entirely repetitive (discarded).\n";
    cerr << nbelowmin << " reads began or were clipped to below the minimum length (discarded).\n";
    cerr << nclipped1side << " reads were clipped on one end.\n";
    cerr << nclipped2side << " reads were clipped on both ends.\n";
    cerr << "Histogram:\n#Bases Clipped\t#Reads\n";
    for (map<int,int>::iterator it=clipped_counts.begin(); it != clipped_counts.end(); ++it) {
        cerr << it->first << "\t" << it->second << endl;
    }
}
