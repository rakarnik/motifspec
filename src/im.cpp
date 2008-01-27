#include "im.h"

int main(int argc, char *argv[]) {
	set_new_handler(out_of_memory);

	string seqfile;                       // file with sequences
	string exprfile;                      // file with expression data
	string clusfile;
	string outfile;
	if(argc < 6) {
    print_usage(cout);
    exit(0);
  }
	
	if ((! GetArg2(argc, argv, "-s", seqfile)) 
				|| (! GetArg2(argc, argv, "-e", exprfile)) 
				|| (! GetArg2(argc, argv, "-o", outfile))) {
		print_usage(cout);
    exit(0);
	}
	
	// Read parameters
	if(! GetArg2(argc, argv, "-numcols", ncol)) ncol = 10;
	if(! GetArg2(argc, argv, "-minsize", minsize)) minsize = 5;
	if(! GetArg2(argc, argv, "-mincorr", mincorr)) mincorr = 0.5;
	
  vector<string> seqs, nameset1;
  cerr << "Reading sequence data from '" << seqfile << "'... ";
	get_fasta_fast(seqfile.c_str(), seqs, nameset1);
	cerr << "done.\n";
	
	vector<string> nameset2;
	cerr << "Reading expression data from '" << exprfile << "'... ";
  expr = get_expr(exprfile.c_str(), &npoints, nameset2);
	cerr << "done.\n";
	
	if(nameset1.size() == nameset2.size()) {
		ngenes = nameset1.size();
	} else {
		cerr << "Inconsistent sizes for sequence and expression datasets\n";
		exit(0);
	}
	
	bool match = true;
	for (int i = 0; i < ngenes; i++) {
		match = match && (nameset1[i] == nameset2[i]);
		if (nameset1[i] != nameset2[i]) {
			cerr << "Row " << i << ": Sequence: '" << nameset1[i] << "'  " << "Expression: '" << nameset2[i] << "'" << endl;
		}
	}
	if(! match) {
		cerr << "Inconsistent gene names for sequence and expression datasets\n";
		exit(0);
	}
	
	cerr << "Successfully read input files -- dataset size is " << ngenes << " genes X " << npoints << " timepoints" << endl;

	
	jcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		jcorr[i] = new float[ngenes];
		for(int j = 0; j < ngenes; j++) {
			jcorr[i][j] = - 2;
		}
	}
	
	jcorr = new float*[ngenes];
	for(int i = 0; i < ngenes; i++) {
		jcorr[i] = new float[ngenes];
		for(int j = 0; j < ngenes; j++) {
			jcorr[i][j] = -2;
		}
	}
	
	cerr << "Reading precomputed correlation values from jcorr.out..." << endl;
	ifstream corrin("jcorr.out");
	int corrcount = 1;
	string line;
	vector<string> values;
	while(getline(corrin, line)) {
		values = split(line, '\t');
		int i = atoi(values[0].c_str());
		int j = atoi(values[1].c_str());
		float jc = atof(values[2].c_str());
		jcorr[i][j] = jcorr[j][i] = jc;
		corrcount++;
		if(corrcount % 1000000 == 0) cerr << "\tRead " << corrcount << " correlation values" << endl;
	}
	cerr << "done." << endl;
	
	cerr << "Setting up SEModel... " << endl;
	SEModel se;
	se.init(seqs, expr, npoints, nameset1, ncol);
	se.modify_params(argc, argv);
	se.set_final_params();
	se.ace_initialize();
	cerr << "done." << endl;
	
	string tmpstr(outfile);
	tmpstr.append(".tmp.ace");
	doit(tmpstr.c_str(), se);
	
	string outstr(outfile);
	outstr.append(".adj.ace");
	ofstream out(outstr.c_str(), ios::trunc);
	print_ace(out, se);
	out.close();
}

void doit(const char* outfile, SEModel& se) {
	for(int g = 0; g < ngenes; g++)
		se.add_possible(g);
	int nruns = 2 * se.possible_positions()
	            / se.get_params().expect
							/ ncol
							/ se.get_params().undersample
							* se.get_params().oversample;
	
	for(int j = 1; j <= nruns; j++) {
		cerr << "\t\tSearch restart #" << j << "/" << nruns << endl;
		se.search_for_motif(minsize, mincorr);
		if(j % 50 == 0) {
			ofstream out(outfile, ios::trunc);
			print_full_ace(out, se);
		}
	}
}

void print_full_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < se.names().size(); x++) out << "#" << x << '\t' << (se.names())[x] << endl;
  out << '\n';
  se.full_output(out);
}

void print_ace(ostream& out, SEModel& se) {
	out << "Parameter values:\n";
  se.output_params(out);
  out << "\nInput sequences:\n";
  for(int x = 0; x < se.names().size(); x++) out << "#" << x << '\t' << se.names()[x] << endl;
  out << endl;
	se.full_output(out);
}

void print_usage(ostream& fout) {
	fout<<"Usage: im -s seqfile -e exprfile (options)\n";
  fout<<" Seqfile must be in FASTA format.\n";
	fout<<" Exprfile must be in tab-delimited format.\n";
  fout<<"Options:\n";
	fout<<" -k          \tnumber of clusters (10)\n";
  fout<<" -numcols    \tnumber of columns to align (10)\n";
  fout<<" -expect     \tnumber of sites expected in model (10)\n";
  fout<<" -gcback     \tbackground fractional GC content of input sequence (0.38)\n";
  fout<<" -minpass    \tminimum number of non-improved passes in phase 1 (200)\n";
  fout<<" -seed       \tset seed for random number generator (time)\n";
  fout<<" -undersample\tpossible sites / (expect * numcols * seedings) (1)\n"; 
  fout<<" -oversample\t1/undersample (1)\n";
}
