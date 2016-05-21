# cuBLASTP
##cuBLASTP:Fine-Grained Parallelization of Protein Sequence Search on CPU+GPU

###Installation

* Download and install Modern GPU library

cuBLASTP uses segmented sort in Modern GPU library by NVlabs. In this version of
cuBLASTP, we uses Modern GPU library version 1.1, which is available for
download at https://github.com/moderngpu/moderngpu/releases/tag/V1.1.

* Build cuBLASTP

Change PATHs in Makefile: Set MGPU_PATH, which refers to the path of Modern GPU
library, and CUDA_PATH, which refers to the location of CUDA toolkit. 

Compile cuBLAST:
```
    $make
```

###Usage 
* Format FASTA database
```
	$ ./formatdb database.fasta
```

* Run cuBLASTP
```
	$ ./cuBLASTP -i query -d database.fasta
```


###Publication
cuBLASTP: Fine-Grained Parallelization of Protein Sequence Search on
CPU+GPU, Jing Zhang, Hao Wang, Wu-chun Feng, Accepted by the IEEE/ACM
Transactions on Computational Biology and Bioinformatics

cuBLASTP: Fine-Grained Parallelization of Protein Sequence Search on a GPU,
Jing Zhang, Hao Wang, Heshan Lin, Wu-chun Feng, In Proceedings of the 28th IEEE
International Parallel and Distributed Processing Symposium (IPDPS), Phoenix, AZ, May
2014.

###License
Please refer to the included LICENSE files in Licensing directory.

###Acknowledgement
This research was supported in part by NSF IIS- 1247693 (BIGDATA) and NSF
CNS-0960081 (MRI), which resulted in the GPU-accelerated HokieSpeed
supercomputer, operated by Advanced Research Computing at Virginia Tech.
