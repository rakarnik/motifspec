#include "standard.h"
#include "seqset.h"
#include "motif.h"

void tseqsetread();
void tmotifconsensus();
void tmotifread();
void tmotifaddcoltoright();
void tmotifaddcoltoleft();

int main(int argc, const char** argv) {
	tseqsetread();
	tmotifconsensus();
	tmotifread();
	tmotifaddcoltoright();
	tmotifaddcoltoleft();
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

void tmotifconsensus() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	Motif m(s, 9, pseudo);
	m.add_site(0, 6, 1);
	m.add_site(1, 9, 1);
	m.add_site(2, 13, 0);
	m.add_site(3, 16, 1);
	m.add_site(4, 8, 0);
	string cons("ACCGTTTCC");
	assert(cons.compare(m.consensus()) == 0);
	cerr << "Passed tmotifconsensus!\n";
}

void tmotifread() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	Motif m(s, 12, pseudo);
	ifstream motin("test.mot");
	m.read(motin);
	motin.close();
	
	assert(m.ncols() == 9);
	string cons("ACCGTTTCC");
	assert(cons.compare(m.consensus()) == 0);
	
	cerr << "Passed tmotifread!\n";
}

void tmotifaddcoltoright() {
	vector<string> seqs;
	vector<string> names;
	get_fasta_fast("test.seq", seqs, names);
	Seqset s(seqs);
	vector<double> pseudo(4);
	pseudo[0] = pseudo[1] = pseudo[2] = pseudo[3] = 0.25;
	Motif m(s, 12, pseudo);
	ifstream motin("test.mot");
	m.read(motin);
	motin.close();

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
	Motif m(s, 12, pseudo);
	ifstream motin("test.mot");
	m.read(motin);
	motin.close();
	
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
	
	cerr << "Passed taddcoltoright!\n";
}