#include "standard.h"
#include "seqset.h"
#include "motif.h"

void tseqsetread();
void tmotifconsensus();
void tmotifread();
void tmotifaddcoltoright();
void tmotifaddcoltoleft();
void tmotifcolumnsample1();
void tmotifcolumnsample2();

int main(int argc, const char** argv) {
	tseqsetread();
	tmotifread();
	tmotifconsensus();
	tmotifread();
	tmotifaddcoltoright();
	tmotifaddcoltoleft();
	tmotifcolumnsample1();
	tmotifcolumnsample2();
	return 0;
}

void tseqsetread() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	for(int i = 0; i < 5; i++)
		assert(s.len_seq(0) == 30);
	cerr << "Passed tseqsetread!\n";
}

void tmotifread() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 12, pseudo, backfreq);
	ifstream motin("test.mot");
	assert(motin.good());
	m.read(motin);
	motin.close();

	assert(m.ncols() == 9);
	assert(m.get_width() == 9);
	string cons("ACCGTTTCC");
	assert(cons.compare(m.consensus()) == 0);
	cerr << "Passed tmotifread!\n";
}

void tmotifconsensus() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 9, pseudo, backfreq);
	m.add_site(0, 6, 1);
	m.add_site(1, 9, 1);
	m.add_site(2, 13, 0);
	m.add_site(3, 16, 1);
	m.add_site(4, 8, 0);
	string cons("ACCGTTTCC");
	assert(cons.compare(m.consensus()) == 0);
	cerr << "Passed tmotifconsensus!\n";
}

void tmotifaddcoltoright() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 9, pseudo, backfreq);
	m.add_site(0, 6, 1);
	m.add_site(1, 9, 1);
	m.add_site(2, 13, 0);
	m.add_site(3, 16, 1);
	m.add_site(4, 8, 0);
	
	// Check sites
	const vector<Site>& sites = m.sites();
	assert(sites[0].posit() == 6);
	assert(sites[1].posit() == 9);
	assert(sites[2].posit() == 13);
	assert(sites[3].posit() == 16);
	assert(sites[4].posit() == 8);

	m.add_col(10);
	
	// Check columns
	for(int i = 0; i < 9; i++)
		assert(m.column(i) == i);
	assert(m.column(9) == 10);
	
	// Check sites
	assert(sites[0].posit() == 6);
	assert(sites[1].posit() == 9);
	assert(sites[2].posit() == 11);
	assert(sites[3].posit() == 16);
	assert(sites[4].posit() == 6);
	
	cerr << "Passed tmotifaddcoltoright!\n";
}

void tmotifaddcoltoleft() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 9, pseudo, backfreq);
	m.add_site(0, 6, 1);
	m.add_site(1, 9, 1);
	m.add_site(2, 13, 0);
	m.add_site(3, 16, 1);
	m.add_site(4, 8, 0);
	
	// Check sites
	const vector<Site>& sites = m.sites();
	assert(sites[0].posit() == 6);
	assert(sites[1].posit() == 9);
	assert(sites[2].posit() == 13);
	assert(sites[3].posit() == 16);
	assert(sites[4].posit() == 8);
	
	m.add_col(-4);
	
	// Check columns
	assert(m.column(0) == 0);
	for(int i = 1; i < 10; i++)
		assert(m.column(i) == i + 3);
	
	// Check sites
	assert(sites[0].posit() == 2);
	assert(sites[1].posit() == 5);
	assert(sites[2].posit() == 13);
	assert(sites[3].posit() == 12);
	assert(sites[4].posit() == 8);
	
	cerr << "Passed tmotifaddcoltoleft!\n";
}

void tmotifcolumnsample1() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 9, pseudo, backfreq);
	m.add_site(0, 6, 1);
	m.add_site(1, 9, 1);
	m.add_site(2, 13, 0);
	m.add_site(3, 16, 1);
	m.add_site(4, 8, 0);	
	assert(m.column_sample() == false);
	cerr << "Passed tmotifcolumnsample1!\n";
}

void tmotifcolumnsample2() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	vector<double> backfreq(4);
	backfreq[0] = backfreq[1] = backfreq[2] = backfreq[3] = 0.25;
	Motif m(s, 9, pseudo, backfreq);
	m.add_site(0, 9, 1);
	m.add_site(1, 12, 1);
	m.add_site(2, 10, 0);
	m.add_site(3, 19, 1);
	m.add_site(4, 5, 0);
	cerr << m.consensus() << '\n';
	m.column_sample();
	m.column_sample();
	const vector<Site>& sites = m.sites();
	assert(sites[0].posit() == 6);
	assert(sites[1].posit() == 9);
	assert(sites[2].posit() == 13);
	assert(sites[3].posit() == 16);
	assert(sites[4].posit() == 8);
	cerr << "Passed tmotifcolumnsample2!\n";
}

