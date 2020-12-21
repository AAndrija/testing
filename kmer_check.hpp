#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <fstream>
#include <map>
#include <set>
#include <tuple>
#include "fasta_parser.hpp"
#include "fastq_parser.hpp"
#include "sequence.hpp"
#include "pink_alignment.hpp"
#include "pink_minimizers.hpp"
#include "Pink_mapperConfig.h"
using namespace std;

//decode kmer
char decode(unsigned int val, bool rev){
	if(val==0) return rev ? 'C' : 'G';
	if(val==1) return rev ? 'A' : 'T';
	if(val==2) return rev ? 'T' : 'A';
	if(val==3) return rev ? 'G' : 'C';
} 

//check presence of kmer in actual genome
void min_check(tuple<unsigned int,unsigned int, bool> kmer, char* ref_list, unsigned int k){
	char chain[k];
	unsigned int maskica = 3;
	unsigned int val;
	for(int i=0;i<k;i++){
		val = 0;
		val = get<0>(kmer) & maskica;
		val = val >> i*2;
		chain[k-1-i] = decode(val, get<2>(kmer));
		maskica = maskica << 2;
	}
	cout << "=========================================\ndecoding kmer: " << get<0>(kmer)
		<< ", " << get<1>(kmer)
		<< ", " << get<2>(kmer)
		<< "\n";
	for(int i=0;i<k;i++){
		cout << chain[i];
	}
	cout << "\nexpected value:\n";
	for(int i=0;i<k;i++){
		cout << ref_list[get<1>(kmer)+i];
	}
	cout << "\n=========================================\n";
	return;
}
