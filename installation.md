# Installation

## Install dependencies:

* Install zlib:
```
	sudo apt-get install libtool autoconf wget perl subversion build-essential 
	sudo apt-get install gfortran libatlas-dev  libatlas-base-dev git
	sudo wget http://www.zlib.net/zlib-1.2.11.tar.gz
	tar zxvf zlib-1.2.11.tar.gz
	cd zlib-1.2.11
	./configure
	make
	sudo make install
```
* Install alsa-lib:
```
	wget ftp://ftp.alsa-project.org/pub/lib/alsa-lib-1.1.4.1.tar.bz2
	tar -jxvf alsa-lib-1.1.4.1.tar
	cd alsa-lib-1.1.4.1 && ./configure
	make && make install
```

## Corss compile Kaldi:
The whole Kaldi is too big for firefly-3399, so cross compiling is recommended. 
### In x86_64 linxu:
* Install cross compile toolchain.
```
	sudo apt-get install gcc-arm-linux-gnueabihf -y
	sudo apt-get install g++-arm-linux-gnueabihf -y
	sudo apt-get install gfortran-arm-linux-gnueabihf -y
	sudo apt-get install gcc-aarch64-linux-gnu binutils-aarch64-linux-gnu  -y
	sudo apt-get -y install gfortran-aarch64-linux-gnu g++-aarch64-linux-gnu
	sudo apt-get install gcc-arm-linux-gnueabi binutils-arm-linux-gnueabi -y
	sudo apt-get install g++-arm-linux-gnueabi gfortran-arm-linux-gnueabi -y
```

* Download Kaldi.
```
	git clone https://github.com/kaldi-asr/kaldi.git
	cd kaldi/tools
```

* Modify Makefile.
```
	-CXX = g++
	-CC = gcc         # used for sph2pipe
	+CXX = aarch64-linux-gnu-g++-5
	+CC = aarch64-linux-gnu-gcc-5          # used for sph2pipe
```
* Cross compile.
```
	make -j4
	cd ../src
	./configure --static --static-fst --openblas-root=../tools/OpenBLAS/install/ --host=aarch64-linux-gnu --use-cuda=no
	make
```
### In firefly-3399:
* Creat directory and copy necessary executable files from x86_64.
```
	export USR_NAME="host_name_in_x86_64"
	export USR_IPADD="IP_address_in_x86_64"
	export KALDI_PATH="absolute_path_of_corss_compiled_Kaldi_in_x86_64"
	sudo mkdir -p kaldi/src/{featbin,bin,gmmbin,latbin}
	sudo scp ${USR_NAME}@${USR_IPADD}:${KALDI_PATH}/src/bin/am-info kaldi/src/bin/
	sudo scp ${USR_NAME}@${USR_IPADD}:${KALDI_PATH}/src/featbin/{add-deltas,apply-cmvn,compute-cmvn-stats,compute-mfcc-feats,copy-feats,modify-cmvn-stats,splice-feats,transform-feats} kaldi/src/featbin/
	sudo scp ${USR_NAME}@${USR_IPADD}:${KALDI_PATH}/src/gmmbin/gmm-latgen-faster kaldi/src/gmmbin
	sudo scp ${USR_NAME}@${USR_IPADD}:${KALDI_PATH}/src/latbin/{lattice-add-penalty,lattice-best-path,lattice-mbr-decode,lattice-prune,lattice-scale} kaldi/src/latbin/
```

## SpeechRecognition Compile
Modify SpeechRecognition/conf/config.hpp, change the KALDI_ROOT to Kaldi root directory in your computer and GRAPH_DIR to model graph director in which the HCLG.fst is. Except the path, there are some other options which have same funtions in Kaldi. To achive a better accuracy, you should adjust the microphone's gain to a proper level.
```
	git clone gitlab@219.139.34.186:openailab/SpeechRecognition.git
	cd SpeechRecognition
	make
	./app
```

## Train
Download data and recipe for the pre-trained model in this project from [CR_resource](ftp://ftp.openailab.net/CR_resource/), or build your own model. After model generated, copy the model files (like, exp/tri1) to command_recognition/app directory. Then modefy the config.hpp, compile and run.

This project only support GMM-HMM acoustic model with input of MFCC feature now. 

Also see http://www.kaldi-asr.org/doc/kaldi_for_dummies.html for data preparing and model training.
