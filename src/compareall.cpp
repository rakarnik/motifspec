#include "standard.h"
#include "archivesites.h"

int main(int argc, char** argv) {
	string seqfile(argv[1]);
	string acefile(argv[2]);
	
	vector<string> seqs;
	get_fasta_fast(seqfile.c_str(), seqs);
	Seqset seqset(seqs);
	BGModel bgm(seqset, 0);
	vector<double> backfreq(4);
	backfreq[0] = backfreq[3] = (1 - bgm.gcgenome())/2.0;
	backfreq[1] = backfreq[2] = bgm.gcgenome()/2.0;
	vector<double> pseudo(backfreq);
	ArchiveSites archive(seqset, bgm, 0.8, pseudo, backfreq);
	ifstream acestream(acefile.c_str());
	archive.read(acestream);
	acestream.close();
	MotifCompare mc(seqset, bgm);
	
	vector<Motif> mots(archive.get_archive());
	int nmots = mots.size();
	for(int i = 0; i < nmots; i++) {
		for(int j = i + 1; j < nmots; j++) {
			cout << acefile << '\t' << i + 1 << '\t' << acefile << '\t' << j + 1 << '\t' << mc.compare(mots[i], mots[j]) << '\n';
		}
	}
}
