#include "conf/config.hpp"
#include <alsa/asoundlib.h>
#include <mutex>

using namespace std;
#define FIFO_LEN 1000000
typedef pair<int, map<int, int>> MAP_PAIR;

typedef struct PCM_params
{
	snd_pcm_format_t SampleFormat = SND_PCM_FORMAT_S16_LE;
	unsigned short nChannleNumber = 2;
	unsigned int nSampleRate = 16000;       //may changed in DeviceInit
	unsigned long nBytesPerSecond = 16000*2*16/8; //SampleRate * NumChannels * BitsPerSample/8, may changed in DeviceInit
	unsigned short BlockAlign = 2 * 16 / 8;  //NumChannels*BitsPerSample/8
	unsigned long nBitsPerSample = 16;       //
	unsigned short nAudioFormat = 1;         //pcm format	
} PCM_params;

typedef struct HXD_WAVFLIEHEAD
{
	char RIFFNAME[4];
	uint32_t nRIFFLength;
	char WAVNAME[4];
	char FMTNAME[4];
	uint32_t nFMTLength;
	uint16_t nAudioFormat;
	uint16_t nChannleNumber;
	uint32_t nSampleRate;
	uint32_t nBytesPerSecond;
	uint16_t BlockAlign;
	uint16_t nBitsPerSample;
	char DATANAME[4];
	uint32_t nDataLength;
} HXD_WAVFLIEHEAD;

class PCM_FIFO
{
public:
	PCM_FIFO() = default;
	
	void write_fifo(const char *data, uint32_t len);
	uint32_t read_fifo(char *data, uint32_t len);
	uint32_t num_in_fifo(void);
private:
	char databuf[FIFO_LEN];
	uint32_t readptr = 0;
	uint32_t writeptr = 0;
	uint32_t num;
};

int CheckExist(const string &filename);
int CreateDir(const string &sPathName);
int KaldiDataPre(const string &path, const string &wav_name);
int MakeMfcc(const string &data_dir, const string &log_dir, const string &mfcc_dir);
int ComputeCmvnStats(const string &data_dir, const string &log_dir, const string &cmvn_dir);
int Decode(const string &graph_dir, const string &data_dir, const string &decode_dir);
int Split(const string& src, const string &separator, vector<string> &dest);
void TerminalCall(const string &command, const string &from);
bool CmpBySumValue(const MAP_PAIR &l, const MAP_PAIR &r);
int MapSum(const map<int, int>& map_to_compute);
float GetMean(const map<int, int> &map_to_compute);
int GetPercentile(const map<int, int> &length_to_count, float fraction);
void ExitWithInfo(const string & info, const string& from);
void OutputInfo(const string & info, const string& from);
void CallBack(const string &result);
int DeviceCapture(snd_pcm_t *pcm_handle, PCM_params &pcm_params, char *data, uint32_t frams_to_record);
bool DeviceInit(const char *pcm_name, snd_pcm_t **pcm_handle, PCM_params *pcm_params);
int CreateWAVHead(const PCM_params &pcm_params, HXD_WAVFLIEHEAD *wav_head, unsigned long len);
bool IsSilence(const char * data, uint32_t len);
void *RecordThread(void);
void *DecideThread(void);
void *DecodeThread(void);

PCM_FIFO pcm_fifo;
static vector<string> record_name;
static PCM_params pcm_params;
mutex fifo_lock, vector_lock;

int main()
{
	thread record_thread(RecordThread);

	thread decide_thread(DecideThread);

	thread decode_thread(DecodeThread);

	record_thread.join();
	decide_thread.join();
	decode_thread.join();
	cout << "error: main thread quit!" << endl;
}


void *RecordThread(void)
{
	snd_pcm_t *pcm_handle;	
	DeviceInit("default", &pcm_handle, &pcm_params);
	uint32_t p_frames = 800;
	char data[4000];   // no less than period frames * BlockAlign
	while(1)
	{
		uint32_t length_readed = DeviceCapture(pcm_handle, pcm_params, data, p_frames);
		pcm_fifo.write_fifo(data, length_readed);
	}
}


void *DecideThread(void)
{
	bool is_silence, need_decode = false, file_opened = false;
	int chunksize = 320;         //chunk length 5ms, 5ms * 64byte/ms = 320bytes
	int delay_threshold = 80;    //only if silence lasted longer than 80 chunks(0.4s), stop to transmit data to decoder
	int delayed_times = 0;
	int filenum = 0, file_length = 0, max_length = 320000; 
	char data[chunksize];
	string filename;
	ofstream wav_out;
	HXD_WAVFLIEHEAD wav_head;
	while(1)
	{
		if(pcm_fifo.num_in_fifo() > chunksize)  
		{
			int real_length = pcm_fifo.read_fifo(data, chunksize);
			is_silence = IsSilence(data, real_length);
			//check whether the data need to decode
			if(!is_silence)
			{
				need_decode = true;
				delayed_times = delay_threshold;
			}
			else if(need_decode == true)
			{
				delayed_times--;
				if(delayed_times <= 0)
				{
					need_decode = false;
				}
			}
			//need transmit and outputfile isn't opened, then open new file
			if(need_decode && !file_opened)
			{
				filenum++;
				filename = string("Test-") + to_string(filenum);
				string realname = string("cache/sound/") + filename + ".wav";
				wav_out.open(realname, ofstream::out | ofstream::ate | ofstream::binary);
				if(!wav_out.is_open())
					ExitWithInfo("wav outpout file open failed", __FUNCTION__);
				CreateWAVHead(pcm_params, &wav_head, 0);        //only for reserve wav-head space in file head
				wav_out.write((char*)(&wav_head), sizeof(HXD_WAVFLIEHEAD));
				file_length = 0;
				file_opened = true;
			}
			//data needed to decode and output file opened, then write data to file
			if(need_decode && file_opened)
			{
				wav_out.write(data, real_length);
				file_length += real_length;
			}
			//need to close the output file and send file to decoder
			if(file_opened && (!need_decode || file_length >= max_length))
			{
				CreateWAVHead(pcm_params, &wav_head, file_length);
				wav_out.seekp(ios_base::beg);
				wav_out.write((char*)(&wav_head), sizeof(HXD_WAVFLIEHEAD));
				wav_out.close();
				file_opened = false;
				vector_lock.lock();
				record_name.push_back(filename);
				vector_lock.unlock();
			}
			
		}
		else
		{
			this_thread::sleep_for(chrono::milliseconds(100));  //sleep 10ms, wait for data
		}
	}
}


void *DecodeThread(void)
{	
/*
	const string filename("Test-0001");
	CreateDir(string("cache/sound"));

	//sounf recording
	TerminalCall(string("rec -r 16000 cache/sound/Test-0001.wav trim 0 00:20"), __FUNCTION__);

	//set dynamic link library path
	char *ld_library = getenv("LD_LIBRARY_PATH");
	if(ld_library != NULL)
		setenv("LD_LIBRARY_PATH", (OPENFST_DIR + ":" + ld_library).c_str(), 1);
	else
		setenv("LD_LIBRARY_PATH", (OPENFST_DIR).c_str(), 1);
	//TerminalCall(string("echo $LD_LIBRARY_PATH"), __FUNCTION__);
*/
	static bool free_last;
	CreateDir("cache/data/test");
	CreateDir("cache/data/mfcc/test");
	CreateDir("cache/exp/make_mfcc/test");
	CreateDir("cache/mfcc/test");
	CreateDir("cache/exp/mfcc_cmvn/test");
	CreateDir("cache/log");
	CreateDir(DECODE_DIR);
	while(1)
	{
		vector_lock.lock();
		int size = record_name.size();
		vector_lock.unlock();
		if(size > 0)
		{
			free_last = false;
			cout << "Recognizing ..." << endl;
			vector_lock.lock();
			string filename = record_name.front();
			record_name.erase(record_name.begin());
			vector_lock.unlock();
			struct timeval tv;
			gettimeofday(&tv,NULL); 
			double start_time = tv.tv_sec * 1000 + tv.tv_usec/1000;

			KaldiDataPre(string("cache/data/test"),filename);	

			KaldiDataPre(string("cache/data/mfcc/test"),filename);

			MakeMfcc("cache/data/mfcc/test", "cache/exp/make_mfcc/test", "cache/mfcc/test");

			ComputeCmvnStats("cache/data/mfcc/test", "cache/exp/mfcc_cmvn/test", "cache/mfcc/test");

			if(Decode(GRAPH_DIR, "cache/data/mfcc/test", DECODE_DIR) == 1)
				cout << "No word match" << endl;

			gettimeofday(&tv,NULL);
			double end_time = tv.tv_sec * 1000 + tv.tv_usec/1000;
			OutputInfo("decode time is : " + to_string(end_time - start_time) + " ms", __FUNCTION__);

			remove((string("cache/sound/") + filename + ".wav").c_str());    //remove decoded wav file
		}
		else
		{
			this_thread::sleep_for(chrono::milliseconds(10));           //sleep 10ms, wait for file
			if(!free_last)
			{
				cout << "Waiting for command" << endl;
				free_last = true;
			}
		}
	}
}


int CreateDir(const string &sPathName)
{
	string DirName(sPathName);
	int i = 0, len = 0;
	
	if(*DirName.end() != '/')
		DirName.push_back('/');
	len = DirName.size();
	if(access(DirName.c_str(),0) == 0)
		return 0;
	for(i = 1; i < len; i++)
	{
		if(DirName[i] == '/')
		{	
			DirName[i] = '\0';
			if(access(DirName.c_str(), 0) != 0)
			{
				if(mkdir(DirName.c_str(), 0777) == -1)  //failed to make dir
					return -1;
			}
			DirName[i] = '/';
		}

	}
	return 0;
}

int KaldiDataPre(const string &path,const string &wav_name)
{
	char pwd[256];
	string spk_name, file_name;

	getcwd(pwd,sizeof(pwd));      //get present work dir

	spk_name = wav_name.substr(0,wav_name.find('-'));

	ofstream utt2ppk_out(path + "/utt2spk");
	utt2ppk_out << wav_name << " " << spk_name <<endl;


	ofstream spk2utt_out(path + "/spk2utt");
	spk2utt_out << spk_name << " " << wav_name << endl;


	ofstream text_out(path + "/text");
	text_out << wav_name << " 0" << endl;


	ofstream wav_out(path + "/wav.scp");
	wav_out << wav_name << " " << pwd << "/cache/sound/" << wav_name << ".wav" << endl;


	ofstream word_out(path + "/words.txt");
	word_out << wav_name << " 0" << endl;
}

int MakeMfcc(const string &data_dir, const string &log_dir, const string &mfcc_dir)
{
	string data_path, scp_path, mfcc_path, segments_path,split_scp, write_num_frames_opt, lines, command;
	char pwd[256];
	getcwd(pwd,sizeof(pwd));      //get present work dir
	data_path = string(pwd) + "/" + data_dir;

	scp_path = data_dir + "/wav.scp";
	if(access(scp_path.c_str(), 0) !=0 )
		ExitWithInfo(scp_path + "not found", __FUNCTION__);

	if(write_utt2num_frames)
		write_num_frames_opt = "--write-num-frames=ark,t:" +log_dir + "/utt2num_frames.1";
	else
		write_num_frames_opt = "\0";
	string spk2warp_path = data_path + "/spk2warp", utt2warp_path = data_path + "/utt2warp", vtln_opts = "\0";
	if (access(spk2warp_path.c_str(),0) == 0)
	{
		OutputInfo(string("using VTLN warp factors from ") + data_dir + "/spk2warp", __FUNCTION__);
		vtln_opts="--vtln-map=ark:" + data_path + "/spk2warp --utt2spk=ark:" + data_path + "/utt2spk";
	}	
	else if (access(utt2warp_path.c_str(),0) == 0)	
	{
		OutputInfo(string("using VTLN warp factors from ") + data_dir + "/utt2warp", __FUNCTION__);
		vtln_opts="--vtln-map=ark:" + data_path + "/utt2warp";
	}

	ifstream wav_in(scp_path);
	if(!wav_in.is_open())
		ExitWithInfo(scp_path + "open failed", __FUNCTION__);
	segments_path = data_dir + "/segments";
	if(access(segments_path.c_str(), 0) == 0)
	{
		ExitWithInfo("segments file exists: not implemented", __FUNCTION__);
	}
	else
	{
		OutputInfo("no segments file exists: assuming wav.scp indexed by utterance.", __FUNCTION__);
		split_scp = log_dir + "/wav_test.1.scp";
		ofstream scp_out(split_scp);
		while(getline(wav_in,lines))
			scp_out << lines << endl;
		command = KALDI_ROOT + "/src/featbin/compute-mfcc-feats " + vtln_opts + "--verbose=2 --config=" \
			+ MFCC_config_path + " scp,p:" + split_scp + " ark:cache/mfcc.ark 2>cache/log/compute-mfcc.log";
/*
		command = KALDI_ROOT + "/src/featbin/compute-mfcc-feats " + vtln_opts + "--verbose=2 --use-energy=" \
			+ use_energy + "--sample-frequency=" + sample_rate  + " scp,p:" + split_scp \
			+ " ark:cache/mfcc.ark 2>cache/log/compute-mfcc.log";
*/
		TerminalCall(command, __FUNCTION__);
		command = KALDI_ROOT + "/src/featbin/copy-feats " +  write_num_frames_opt \
			+ " --compress=" + mfcc_compress + " ark:cache/mfcc.ark ark,scp:" + mfcc_dir \
			+ "/raw_mfcc_test.1.ark," + mfcc_dir + "/raw_mfcc_test.1.scp 2>cache/log/copy-feats.log";
		TerminalCall(command, __FUNCTION__);
	}
	wav_in.close();

	//if error occured
	string error_log = log_dir + "/.error.test";
	if(access(error_log.c_str(), 0) == 0)
		ExitWithInfo(string("error producing mfcc features for test: ") + log_dir \
			     + "/make_mfcc_test.1.log", __FUNCTION__);

	//concatenate the .scp files together, only 1 job
	ofstream feat_out(data_path + "/feats.scp");
	for(int j = 1; j <= nj; j++)
	{
		ifstream scp_in(mfcc_dir + "/raw_mfcc_test." +to_string(j) + ".scp");
		if(!scp_in.is_open())
			ExitWithInfo(mfcc_dir + "/raw_mfcc_test." +to_string(j) + ".scp open failed", __FUNCTION__);
		while(getline(scp_in, lines))
			 feat_out << lines << endl;
	}
	feat_out.close();

	if(write_utt2num_frames)
	{
		ofstream utt2num_out(log_dir + "/utt2num_frames");
		for(int j = 1; j <= nj; j++)
		{
			ifstream utt2num_in(log_dir + "/utt2num_frames." + to_string(j));
			if(!utt2num_in.is_open())
				ExitWithInfo(log_dir + "/utt2num_frames." + to_string(j) + " open failed", __FUNCTION__);
			while(getline(utt2num_in, lines))
				 utt2num_out << lines << endl;
			remove((log_dir + "/utt2num_frames." +to_string(j)).c_str());
		}
	}
	
	int utts_num = 0, feats_num = 0;
	ifstream utt2spk_in(data_path + "/utt2spk");
	if(!utt2spk_in.is_open())
		ExitWithInfo(data_path + "/utt2spk open failed", __FUNCTION__);
	ifstream feats_in(data_path + "/feats.scp");
	if(!feats_in.is_open())
		ExitWithInfo(data_path + "/feats.scp open failed", __FUNCTION__);
	while(getline(utt2spk_in, lines))
		utts_num++;
	while(getline(feats_in, lines))
		feats_num++;
	if(utts_num != feats_num)
	{
		OutputInfo(string("It seems not all of the feature files were successfully processed ") \
			   + to_string(feats_num) + " feats, " + to_string(utts_num) +" utts.", __FUNCTION__);
		if(feats_num < (utts_num - utts_num/20))
			ExitWithInfo("Less than 95% the features were successfully generated.  Probably a serious error."\
				      , __FUNCTION__);
	}
}

int ComputeCmvnStats(const string &data_dir, const string &log_dir, const string &cmvn_dir)
{
 	char pwd[256];
	getcwd(pwd,sizeof(pwd));      //get present work dir
	string cmvn_path, lines, command;
	cmvn_path = string(pwd) + "/" + cmvn_dir;
	
	//check files
	if(access(string(data_dir + "/feats.scp").c_str(), 0) != 0)
		ExitWithInfo(data_dir + "/feats.scp not found", __FUNCTION__);
	if(access(string(data_dir + "/spk2utt").c_str(), 0) != 0)
		ExitWithInfo(data_dir + "/spk2utt not found", __FUNCTION__);
	
	if(cmvn_fake)
	{
		ExitWithInfo("cmvn_fake: not implemented", __FUNCTION__);
	}
	else if(cmvn_two_channel)
	{
		command = KALDI_ROOT + "/src/featbin/compute-cmvn-stats-two-channel " + data_dir \
			+ "/reco2file_and_channel scp:" + "/feats.scp ark,scp:" + cmvn_dir + "/cmvn_test.ark," \
			+ cmvn_dir + "/cmvn_test.scp 2>cache/log/compute_cmvn_states.log"; 
		TerminalCall(command, __FUNCTION__);	
	}
	else if(fake_dims != "\0")
	{
		command = KALDI_ROOT + "/src/featbin/compute-cmvn-stats --spk2utt=ark:" + data_dir \
			+ "/spk2utt scp:" + data_dir + "/feats.scp ark:cache/fake_dims.ark 2>cache/log/compute_cmvn_states.log";
		TerminalCall(command, __FUNCTION__);
		command = KALDI_ROOT + "/src/featbin/modify-cmvn-stats" + fake_dims + "ark:cache/fake_dims.ark ark,scp:" \
			+ cmvn_dir + "/cmvn_test.ark," + cmvn_dir + "/cmvn_test.scp 2>cache/log/modify_cmvn_stats.log";
		TerminalCall(command, __FUNCTION__);
	}
	else
	{
		command = KALDI_ROOT + "/src/featbin/compute-cmvn-stats --spk2utt=ark:" + data_dir \
			+ "/spk2utt scp:" + data_dir + "/feats.scp ark,scp:" + cmvn_dir + "/cmvn_test.ark," \
			+ cmvn_dir + "/cmvn_test.scp 2>cache/log/compute_cmvn_stats.log";	
		TerminalCall(command, __FUNCTION__);
	}
	

	ifstream cmvn_test_in(cmvn_dir + "/cmvn_test.scp");
	if(!cmvn_test_in.is_open())
		ExitWithInfo(cmvn_dir + "/cmvn_test.scp open failed", __FUNCTION__);
	ofstream cmvn_out(data_dir + "/cmvn.scp");
	while(getline(cmvn_test_in, lines))
		cmvn_out << lines << endl;
	cmvn_test_in.close();
	cmvn_out.close();

	int cmvn_num = 0, spk_num = 0;
	ifstream cmvn_in(data_dir + "/cmvn.scp");
	if(!cmvn_in.is_open())
		ExitWithInfo(data_dir + "/cmvn.scp open failed", __FUNCTION__);
	ifstream spk2utt_in(data_dir + "/spk2utt");
	if(!spk2utt_in.is_open())
		ExitWithInfo(data_dir + "/spk2utt open failed", __FUNCTION__);
	while(getline(cmvn_in, lines))
		cmvn_num++;
	while(getline(spk2utt_in, lines))
		spk_num++;	
	if(cmvn_num != spk_num)
		ExitWithInfo(string("compute cmvn: warning: it seems not all of the speakers got cmvn stats ") \
			     + to_string(cmvn_num) + " cmvn, " + to_string(spk_num) + " spk", __FUNCTION__);
	OutputInfo("Succeeded creating CMVN stats for test", __FUNCTION__);
}


//<graph-dir> <data-dir> <decode-dir>
int Decode(const string &graph_dir, const string &data_dir, const string &decode_dir)
{
	string lines, decode_src, feat_type, feats, thread_string, command, features_rspecifier, words;
	vector<string> split_words, results;
	//only one utt, no need to split, just copy files
	CreateDir(data_dir + "/split1/1");
	vector<string> files_copy = {"/cmvn.scp", "/feats.scp", "/spk2utt", "/text", "/utt2spk", "/wav.scp"};
	for(int i =0; i < files_copy.size(); i++)
	{
		ifstream file_in(data_dir + files_copy[i]);
		if(!file_in.is_open())
			ExitWithInfo(data_dir + files_copy[i] + " open failed", __FUNCTION__);
		ofstream file_out(data_dir + "/split1/1" + files_copy[i]);
		if(!file_out.is_open())
			ExitWithInfo(data_dir + "/split1/1" + files_copy[i] + " open failed", __FUNCTION__);
		while(getline(file_in,lines))
			file_out << lines << endl;
	}

	//record number of jobs
	ofstream num_jobs_out(decode_dir + "/num_jobs");
	num_jobs_out << nj << endl;

	decode_src = decode_dir.substr(0, decode_dir.rfind('/'));

	if(decode_model == string("\0"))
	{
		if(decode_iter == string("\0"))
			decode_model = decode_src + "/final.mdl";
		else
			decode_model = decode_src + "/" + decode_iter + ".mdl";
	}
	
	vector<string> files_check = {data_dir + "/split1/1/feats.scp", data_dir + "/split1/1/cmvn.scp", \
				      decode_model, graph_dir + "/HCLG.fst"};
	for(int i = 0; i < files_check.size(); i++)
	{
		if(access(files_check[i].c_str(), 0) != 0)
			ExitWithInfo(files_check[i] + " not found", __FUNCTION__);
	}
	
	if(access((decode_src + "/final.mat").c_str(), 0) == 0)
		feat_type = "lda";
	else
		feat_type = "delta";
	OutputInfo(string("feature type is ") + feat_type, __FUNCTION__);
	
	if(decode_num_threads > 1)
		thread_string = "-parallel --num-threads=$num_threads";
	else
		thread_string = "\0";

	if(decode_stage <= 0)
	{
		if(access((graph_dir + "/num_pdfs").c_str(), 0) ==0)
		{
			command = KALDI_ROOT + "/src/bin/am-info --print-args=false " + decode_model \
				+ " | grep pdfs | awk '{print $NF}' > cache/pdfs_in_model";
			TerminalCall(command, __FUNCTION__);
			ifstream model_pdfs_in("cache/pdfs_in_model");
			if(!model_pdfs_in.is_open())
				ExitWithInfo("cache/pdfs_in_model open failed", __FUNCTION__);
			getline(model_pdfs_in, lines);
			int pdfs_in_model = stod(lines);
			ifstream num_dpfs_in(graph_dir + "/num_pdfs");
			if(!num_dpfs_in.is_open())
				ExitWithInfo(graph_dir + "/num_pdfs open failed", __FUNCTION__);
			getline(num_dpfs_in, lines);
			int pdfs_in_file = stod(lines);
			if(pdfs_in_file != pdfs_in_model)
				ExitWithInfo(string("mismatch in number of pdfs with ") + decode_model, __FUNCTION__);
		}
		
		string cmvn_opts;
		ifstream cmvn_opts_in(decode_src + "/cmvn_opts");
		if(cmvn_opts_in.is_open())
			getline(cmvn_opts_in, cmvn_opts);
		command = KALDI_ROOT + "/src/featbin/apply-cmvn " + cmvn_opts +" --utt2spk=ark:" + data_dir \
			+ "/split1/1/utt2spk scp:" + data_dir + "/split1/1/cmvn.scp scp:" + data_dir \
			+ "/split1/1/feats.scp ark:cache/apply-cmvn.ark 2>cache/log/apply_cmvn.log";
		TerminalCall(command, __FUNCTION__);
			
		if(feat_type == "lda")
		{
			string splice_opts;
			ifstream splice_opts_in(decode_src + "/splice_opts");
			if(splice_opts_in.is_open())
				getline(splice_opts_in, splice_opts);
			command = KALDI_ROOT + "/src/featbin/splice-feats " + splice_opts \
				+ " ark:cache/apply-cmvn.ark" + " ark:cache/splice-feats.ark 2>cache/log/splice_feats.log";
			TerminalCall(command, __FUNCTION__);
			command =  KALDI_ROOT + "/src/featbin/transform-feats " + decode_src \
				+ "/final.mat ark:cache/splice-feats.ark ark:cache/trasform-feats.ark 2>cache/log/transform_feats.log"; 
			TerminalCall(command, __FUNCTION__);
			features_rspecifier = "cache/trasform-feats.ark";
		}
		else if(feat_type == "delta")
		{
			string delta_opts;
			ifstream delta_optts_in(decode_src + "/delta_opts");
			if(delta_optts_in.is_open())
				getline(delta_optts_in, delta_opts);
			command = KALDI_ROOT + "/src/featbin/add-deltas " + delta_opts \
				+ " ark:cache/apply-cmvn.ark ark:cache/add-delta.ark 2>cache/log/add_deltas.log";
			TerminalCall(command, __FUNCTION__);
			features_rspecifier = "cache/add-delta.ark";
		}

		if(transform_dir != "\0")
		{
			OutputInfo(string("using fMLLR transforms from ") + transform_dir, __FUNCTION__);
			if(access((transform_dir + "/trans.1").c_str(), 0) != 0)
				ExitWithInfo(string("Expected ") + transform_dir + "/trans.1 to exist.", __FUNCTION__);
			ifstream trans_num_in(transform_dir + "/num_jobs");
			if(!trans_num_in.is_open())
				ExitWithInfo(transform_dir + "/num_jobs open failed", __FUNCTION__);
			getline(trans_num_in, lines);
			int nj_orig = stod(lines);
			if(nj != nj_orig)
			{
				ofstream all_trans_out("cache/all_tans.ark");
				for(int i = 0; i < nj_orig; i++)
				{
					ifstream trans_in(transform_dir + "/trans." + to_string(i));
					if(!trans_in.is_open())
						ExitWithInfo(transform_dir + "/trans." + to_string(i) + " open failed", __FUNCTION__);
					while(getline(trans_in, lines))
						all_trans_out << lines << endl;
				}
				command = KALDI_ROOT + "/src/featbin/copy-feats ark:cache/all_tans.ark ark,scp:" \
					+ decode_dir + "/trans.ark," + decode_dir + "/trans.scp";
				TerminalCall(command, __FUNCTION__);
				command = KALDI_ROOT + "/src/featbin/transform-feats --utt2spk=ark:" + data_dir \
					+ "/split1/1/utt2spk scp:" + decode_dir + "/trans.scp ark:" + features_rspecifier\
					+ " ark:cache/trasform-feats2.ark 2>cache/log/transform_feats.log";
				TerminalCall(command, __FUNCTION__);
				features_rspecifier = "cache/trasform-feats2.ark";
			}
			else
			{
				command = KALDI_ROOT + "/src/featbin/transform-feats --utt2spk=ark:" + data_dir \
					+ "/split1/1/utt2spk ark:" + decode_dir + "/trans.1.ark ark:" \
					+ features_rspecifier + " ark:cache/trasform-feats2.ark 2>cache/log/transform_feats.log";
				TerminalCall(command, __FUNCTION__);
				features_rspecifier = "cache/trasform-feats2.ark";
			}
		}

		command = KALDI_ROOT + "/src/gmmbin/gmm-latgen-faster" + thread_string + " --max-active=" \
			+ decode_max_active + " --beam=" + decode_beam + " --lattice-beam=" + decode_lattice_beam \
			+ " --acoustic-scale=" + decode_acwt + " --allow-partial=true --word-symbol-table=" + graph_dir \
			+ "/words.txt " + decode_model + " " + graph_dir + "/HCLG.fst ark:" + features_rspecifier \
			+ " ark:- 1>" + decode_dir + "/lat.1 2>cache/log/gmm_latgen_faster.log";
		TerminalCall(command, __FUNCTION__);
		//command = string("gzip -c > ") + decode_dir + "/lat.1.gz";
		//ret = system(command.c_str());
	}

	if(!skip_scoring)
	{
		CreateDir(decode_dir + "/scoring_kaldi/penalty_" + word_ins_penalty);
		vector<string> score_check = {string(graph_dir+"/words.txt"), string(decode_dir + "/lat.1"), string(data_dir + "/text")};
		for(int i; i < score_check.size(); i++)
		{
			if(access(score_check[i].c_str(), 0) != 0)
				ExitWithInfo(score_check[i] + " not found", __FUNCTION__);
		}
		
		{
			ifstream text_in(data_dir + "/text");
			if(!text_in.is_open())
				ExitWithInfo(data_dir + "/text open failed", __FUNCTION__);
			ofstream text_filt(decode_dir + "/scoring_kaldi/test_filt.txt");
			while(getline(text_in,lines))
				text_filt << lines << endl;
		}
		

		if(score_stage <= 0)
		{

			command = KALDI_ROOT + "/src/latbin/lattice-scale --inv-acoustic-scale=" + score_lmwt +" ark:" \
				+ decode_dir +"/lat.1 ark:cache/scaled.lats 2>cache/log/lattice_scale.log";
			TerminalCall(command, __FUNCTION__);

			command = KALDI_ROOT + "/src/latbin/lattice-add-penalty --word-ins-penalty=" + word_ins_penalty \
				+ " ark:cache/scaled.lats ark:cache/add.lats 2>cache/log/lattice_add_penalty.log";
			TerminalCall(command, __FUNCTION__);

			command = KALDI_ROOT + "/src/latbin/lattice-prune --beam=" + score_beam \
				+ " ark:cache/add.lats ark:cache/prune.lats 2>cache/log/latiice_prune.log";
			TerminalCall(command, __FUNCTION__);
			if(score_decode_mbr)
			{
				OutputInfo(string("scoring with MBR, word insertion penalty=") + word_ins_penalty, __FUNCTION__);	
				command = KALDI_ROOT + "/src/latbin/lattice-mbr-decode --word-symbol-table=" + graph_dir \
					+ "/words.txt ark:cache/prune.lats ark,t:cache/final.lats 2>cache/log/lattice_mbr_decode.log";
				TerminalCall(command, __FUNCTION__);
			}
			else
			{
				OutputInfo(string("scoring with word insertion penalty=") + word_ins_penalty, __FUNCTION__);
				command = KALDI_ROOT + "/src/latbin/lattice-best-path --word-symbol-table=" + graph_dir \
					+ "/words.txt ark:cache/prune.lats ark,t:cache/final.lats 2>cache/log/lattice_best_path.log";
				TerminalCall(command, __FUNCTION__);
			}
			map<int, string> int2sym;
			ifstream word_in(graph_dir + "/words.txt");
			if(!word_in.is_open())
				ExitWithInfo(graph_dir + "/words.txt open failed", __FUNCTION__);
			while(getline(word_in,lines))
			{
				split_words.clear();
				Split(lines, " ", split_words);
				if(split_words.size() != 2)
					ExitWithInfo(string("bad line in symbol table file: ") + graph_dir + "/words.txt", __FUNCTION__);
				int2sym.insert(make_pair(stod(split_words[1]), split_words[0]));
			}

			string filename = decode_dir + "/scoring_kaldi/penalty_" + word_ins_penalty + "/" \
					+ score_lmwt + ".txt";
			ofstream LMWT_out(filename);
			ifstream final_lats("cache/final.lats");
			if(!final_lats.is_open())
				ExitWithInfo("cache/final.lats open failed", __FUNCTION__);
			while(getline(final_lats, lines))
			{
				split_words.clear();
				Split(lines, " ", split_words);
				for(int i = 0; i < split_words.size(); i++)
				{
					if((score_field_begin == "\0" || i >= stod(score_field_begin))&&\
					   (score_field_end == "\0" || i <= stod(score_field_end)))
					{
						int num = stod(split_words[i]);
						auto iter = int2sym.find(num);
						if(iter->second.substr(0, 1) != "<")    //not "<XXX>" style
							results.push_back(iter->second);
					}
				}
				words = "\0";
				for(int i=0; i<results.size(); i++)
				{
					words = words + results[i] + " ";
				}
				if(words != "\0")
				{
					CallBack(words);
					return 0;  //decode success and words match
				}
				else
					return 1; //no word match
			}
		}		
	}
	return -1;//decode failed
}

int Split(const string& src, const string& separator, vector<string>& dest)
{
    string str = src;
    string substring;
    string::size_type start = 0, index;
    start = str.find_first_not_of(separator,start);
    do
    {
        index = str.find_first_of(separator,start);
        if (index != string::npos)
        {    
            substring = str.substr(start,index-start);
            dest.push_back(substring);
            start = str.find_first_not_of(separator,index);
            if (start == string::npos) 
		return -1;
        }
    }while(index != string::npos);
    
    //the last token
    substring = str.substr(start);
    dest.push_back(substring);
}

void TerminalCall(const string &command, const string &from)
{
	OutputInfo(command, from);
	int ret = system(command.c_str());
	if(ret !=0 )
		ExitWithInfo(command + " call in terminal failed.", from);
}

int GetPercentile(const map<int, int> &length_to_count, float fraction)
{
	int total_phones = 0;
	for(auto i:length_to_count)
	{
		total_phones += i.second;
	}
	if(total_phones == 0)
		return 0;
	else
	{
		int count_cutoff = fraction * total_phones;
		int cur_count_total = 0;
		for(auto j:length_to_count)
		{
			//if(j.second <= 0)
			//exit   
			cur_count_total += j.second;
			if (cur_count_total >= count_cutoff)
				return j.first;
		}
	}
}

bool CmpBySumValue(const MAP_PAIR &l, const MAP_PAIR &r)
{
	int l_sum = MapSum(l.second);
	int r_sum = MapSum(r.second);
	return l_sum > r_sum;
}

int MapSum(const map<int, int> &map_to_compute)
{
	int sum = 0;
	for(auto i:map_to_compute)
	{
		sum += i.second;
	}
	return sum;
}

float GetMean(const map<int, int> &map_to_compute)
{
	int total_phones = MapSum(map_to_compute);
	if(total_phones == 0)
		return 0;
	float total_frames = 0;
	for(auto i:map_to_compute)
	{
		total_frames += i.first * i.second;
	}
	return total_frames / total_phones;
}

void ExitWithInfo(const string & info, const string& from)
{
	cout << "<" << from << ">: " << info << endl;
	exit(-1);
}

void OutputInfo(const string & info, const string& from)
{
	static ofstream info_out("cache/log/Output.log");
	info_out << "<" << from << ">: " << info << endl;
}

int DeviceCapture(snd_pcm_t *pcm_handle, PCM_params &pcm_params, char *data, uint32_t frams_to_record)
{ 
	snd_pcm_sframes_t len_recorded = 0;

	len_recorded = snd_pcm_readi(pcm_handle, data, frams_to_record);

	if (len_recorded == -EPIPE) 
	{
		/* EPIPE means overrun */
		OutputInfo("alsa overrun occurred", __FUNCTION__);
		snd_pcm_prepare(pcm_handle);
	} 
	else if (len_recorded < 0)
	{
		OutputInfo(string("error: alsa read ") + snd_strerror(len_recorded), __FUNCTION__);
	}
	return len_recorded * pcm_params.BlockAlign;
} 


bool DeviceInit(const char *pcm_name, snd_pcm_t **pcm_handle, PCM_params *pcm_params) 
{
	int err = 0, val = 0;
	//snd_pcm_uframes_t frames;
	 snd_pcm_hw_params_t *hw_params;

	//open device
	if ((err = snd_pcm_open(pcm_handle, pcm_name, SND_PCM_STREAM_CAPTURE, 0)) < 0) 
	{
		ExitWithInfo(string("cannot open audio device") + pcm_name + snd_strerror(err), __FUNCTION__);
	}

	//malooc params space
	if ((err = snd_pcm_hw_params_malloc(&hw_params)) < 0) 
	{
		ExitWithInfo(string("cannot allocate hardware parameter structure ") + snd_strerror(err), __FUNCTION__);
	}

	//init parameters
	if ((err = snd_pcm_hw_params_any(*pcm_handle, hw_params)) < 0) 
	{
		ExitWithInfo(string("cannot initialize hardware parameter structure ") + snd_strerror(err), __FUNCTION__);
	}

	//set access mode
	if ((err = snd_pcm_hw_params_set_access(*pcm_handle, hw_params,\
						SND_PCM_ACCESS_RW_INTERLEAVED)) < 0) 
	{
		ExitWithInfo(string("cannot set access type ") + snd_strerror (err), __FUNCTION__);
	}

	//set PCM sample format 
	if ((err = snd_pcm_hw_params_set_format(*pcm_handle, hw_params, pcm_params->SampleFormat)) < 0)
	{
		ExitWithInfo(string("cannot set sample format ") + snd_strerror (err), __FUNCTION__);
	}

	//set sample rate
	if ((err = snd_pcm_hw_params_set_rate_near(*pcm_handle, hw_params, &pcm_params->nSampleRate, &val)) < 0)
	{
		ExitWithInfo(string("cannot set sample rate ") + snd_strerror (err), __FUNCTION__);
	} 

	//set channels  (two channel)
	if ((err = snd_pcm_hw_params_set_channels(*pcm_handle, hw_params, pcm_params->nChannleNumber)) < 0) 
	{
		ExitWithInfo(string("cannot set channel count ") + snd_strerror (err), __FUNCTION__);
	}
 
	//set parameters
	if ((err = snd_pcm_hw_params(*pcm_handle, hw_params)) < 0) 
	{
		ExitWithInfo(string("cannot set parameters ") + snd_strerror (err), __FUNCTION__);
	}

	//calculate nBytesPerSecond once agine
	pcm_params->nBytesPerSecond = pcm_params->nSampleRate * pcm_params->nChannleNumber \
				    * pcm_params->nBitsPerSample / 8;   //SampleRate * NumChannels * BitsPerSample/8

	//free the parameters 
	snd_pcm_hw_params_free(hw_params);
	return true; 
} 

int CreateWAVHead(const PCM_params &pcm_params, HXD_WAVFLIEHEAD *wav_head, unsigned long len)
{
	wav_head->RIFFNAME[0] = 'R';
	wav_head->RIFFNAME[1] = 'I';
	wav_head->RIFFNAME[2] = 'F';
	wav_head->RIFFNAME[3] = 'F';
	wav_head->nRIFFLength = 36 + len;  //  length of total file (without RIFFNAME and nRIFFLength)

	wav_head->WAVNAME[0] = 'W';
	wav_head->WAVNAME[1] = 'A';
	wav_head->WAVNAME[2] = 'V';
	wav_head->WAVNAME[3] = 'E';

	wav_head->FMTNAME[0] = 'f';
	wav_head->FMTNAME[1] = 'm';
	wav_head->FMTNAME[2] = 't';
	wav_head->FMTNAME[3] = ' ';
	
	wav_head->DATANAME[0] = 'd';
	wav_head->DATANAME[1] = 'a';
	wav_head->DATANAME[2] = 't';
	wav_head->DATANAME[3] = 'a';
	wav_head->nFMTLength = 16;

	wav_head->nAudioFormat    = pcm_params.nAudioFormat;
	wav_head->nChannleNumber  = pcm_params.nChannleNumber;
	wav_head->nSampleRate     = pcm_params.nSampleRate;
	wav_head->nBytesPerSecond = pcm_params.nBytesPerSecond;
	wav_head->BlockAlign      = pcm_params.BlockAlign;
	wav_head->nBitsPerSample  = pcm_params.nBitsPerSample;
	wav_head->nDataLength     = len;
}


void PCM_FIFO::write_fifo(const char *data, uint32_t len)
{
	fifo_lock.lock();
	uint32_t fifo_len = this->num;
	fifo_lock.unlock();
	if((fifo_len + len) > FIFO_LEN){
		OutputInfo("FIFO over flow", __FUNCTION__);
                return;
        }
	if((len + this->writeptr) <= FIFO_LEN)
		memcpy(&this->databuf[this->writeptr], data, len);
	else
	{
		unsigned int back = FIFO_LEN - writeptr, head = len - back;
		memcpy(&this->databuf[this->writeptr], data, back);
		memcpy(&this->databuf[0], data + back, head);
	}
	this->writeptr = (this->writeptr + len) % FIFO_LEN;
	fifo_lock.lock();
	this->num += len; 
	fifo_lock.unlock();
}

uint32_t PCM_FIFO::read_fifo(char *data, uint32_t len)
{
	fifo_lock.lock();
	if(len > this->num)    //if no enough data in fifo
		len = this->num;
	fifo_lock.unlock();
	if((len + this->readptr) <= FIFO_LEN)
		memcpy(data, &this->databuf[this->readptr], len);
	else
	{
		unsigned int back = FIFO_LEN - readptr, head = len - back;
		memcpy(&this->databuf[this->writeptr], data, back);
		memcpy(&this->databuf[0], data + back, head);
	}	
	this->readptr = (this->readptr + len) % FIFO_LEN;
	fifo_lock.lock();
	this->num -= len;
	fifo_lock.unlock();
	return len;
}

uint32_t PCM_FIFO::num_in_fifo(void)
{
	return this->num;
}

//check is the chunk silence, use treshhold of mean absolute value of PCM 
//only for SND_PCM_FORMAT_S16_LE format
bool IsSilence(const char * data, uint32_t len)
{
	short threshold = 0x0300;  //need to find a better value!!
	uint64_t sum = 0;
	short value;
	for(uint32_t i=0; i<(len/2); i++)
	{
		value = ((short*)data)[i];
		if(value < 0)
			value *= -1;
		sum += value;
	}
	if((2*sum/len) > threshold)  //if mean absolute value of PCM biger than threshold
		return false;
	else
		return true;
}

void CallBack(const string &result)
{
        //Call back function, you can modify it 
        cout << "The result is " << result << endl;
 
}

