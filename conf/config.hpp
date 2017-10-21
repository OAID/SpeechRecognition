#ifndef _CONFIG_HPP_
#define _CONFIG_HPP_

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <assert.h>
#include <sys/time.h>
#include <thread>

//#define ANALYZE_LATS 0

//path configurations
std::string KALDI_ROOT  = "/home/firefly/kaldi";
std::string GRAPH_DIR   = "tri3b/graph_word";
std::string DECODE_DIR  = "tri3b/decode_test_word";
std::string OPENFST_DIR = KALDI_ROOT + "/tools/openfst/lib";
				

int nj = 1;

//MFCC configurations
std::string MFCC_config_path  = "conf/mfcc.conf";
bool write_utt2num_frames     = false;
std::string mfcc_compress     = "true";
std::string use_energy        = "false";
std::string sample_rate       = "16000";

//CMVN configurations
bool cmvn_fake          = false;
std::string fake_dims   = "\0";
bool cmvn_two_channel   = false;
//split configurations
bool split_per_spk      = true; 
std::string per_utt     = "\0";                     //"--per-utt";
 
//decode configurations
std::string transform_dir="\0";                 // this option won't normally be used, but it can be used if you want to
                                                // supply existing fMLLR transforms when decoding.
std::string decode_iter  = "\0";
std::string decode_model ="\0";                 // You can specify the model to use
int decode_stage                 = 0;
std::string decode_max_active    = "7000";
std::string decode_beam          = "13.0";
std::string decode_lattice_beam  = "6.0";
std::string decode_acwt          = "0.083333";           // note: only really affects pruning (scoring is on lattices).
int decode_num_threads           = 1;        		// if >1, will use gmm-latgen-faster-parallel
std::string decode_parallel_opts = "\0";        // ignored now.
std::string decode_scoring_opts  = "\0";
bool skip_scoring                = false;
std::string splice_opts          = "\0";
std::string cmvn_opts            = "\0";
std::string delta_opts           = "\0";


//diagnose configuration
float frequency_cutoff_percentage = 0.5;

//score configuration
int  score_stage             = 0;
bool score_decode_mbr        = false;
bool score_stats             = true;
std::string score_beam       = "6"; 			
std::string word_ins_penalty = "1.0";		//please set this same as your best wer model
std::string score_lmwt       = "17";            //

std::string score_field_begin = "1";            //range of score field, don't modify
std::string score_field_end   = "\0";		//range of score field, don't modify
#endif
