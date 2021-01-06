# MMpred
### Distance-assisted multimodal conformation sampling for de novo protein structure prediction



**Developer:**   
                Kailong Zhao
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
		Email: zhaokailong@zjut.edu.cn  

**Contact:**  
                Guijun Zhang, Prof  
                College of Information Engineering  
                Zhejiang University of Technology, Hangzhou 310023, China  
                Email: zgj@zjut.edu.cn  

## 1. INSTALLATION
Binaries for Linux 64 bit system has been included in the package. The Linux binary was compiled using GCC 5.4.0. Users need to have these versions of GCC compilers when using binaries.

Please Follow the below steps to install and configure MMpred:

- Download Rosetta3.10 source package from https://www.rosettacommons.org/software/.
(where `$ROSETTA3`=path-to-Rosetta)

- Copy and paste ``"ClassicAbinitio.cc"`` and ``"ClassicAbinitio.hh"`` from ``"code/"`` folder in MMpred package to ``"$ROSETTA3/main/source/src/protocols/abinitio/"`` folder in Rosetta.

- Compile MMpred source code using the following commands:

```
 $> cd $ROSETTA3/main/source/
 $> ./scons.py AbinitioRelax -j<NumOfJobs> mode=release bin
```

## 2. INPUT
MMpred requires four files to generate models:

	fasta				: fasta file
	distance			: distance map file
	aat000_03_05.200_v1_3		: fragment library with fragment lenth 3
	aat000_09_05.200_v1_3		: fragment library with fragment lenth 9

## 3. RUN
Please follow the below steps to run MMpred:

- Go to the ``"example/"`` folder of MMpred.

- Run MMpred with the following command:

```
 $> $ROSETTA3/main/source/bin/AbinitioRelax.default.linuxgccrelease @flags
```

- Five models are generated in the ``"output_files/"`` folder.


## 4. OUTPUT
Output files of MMpred are stored in the ``"example/output_files/"`` folder, including five predicted models (modelX.pdb).

	model1.pdb
	model2.pdb
	model3.pdb
	model4.pdb
	model5.pdb


## 5. DISCLAIMER
The executable software and the source code of MMpred is distributed free of charge 
as it is to any non-commercial users. The authors hold no liabilities to the performance 
of the program.
