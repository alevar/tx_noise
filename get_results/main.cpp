#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "arg_parse.h"
#include "QuantRes.h"

enum Opt {READS_FP   = 'i',
        STRG_FP       = 's',
        SLMN_FP       = 'm',
        KLST_FP    = 'k',
        OUTPUT_FP      = 'o',
        ANNOTATION_FP  = 'a',
        REAL_BASE = 'r',
        SPLICING_BASE = 'p',
        INTRONIC_BASE = 'n',
        INTERGENIC_BASE = 'g'};

int main(int argc, char** argv) {

    ArgParse args("gtex_cor");
    args.add_string(Opt::READS_FP,"reads","","File containing real all simulated with polyester",true);
    args.add_string(Opt::STRG_FP,"stringtie","","File containing final results of running stringtie",true);
    args.add_string(Opt::SLMN_FP,"salmon","","File containing final results of running salmon",true);
    args.add_string(Opt::KLST_FP,"kallisto","","File containing final results of running kallisto",true);
    args.add_string(Opt::ANNOTATION_FP,"gff","","annotation corresponding to the simulated reads",true);
    args.add_string(Opt::REAL_BASE,"real","","base name for real subset",false);
    args.add_string(Opt::SPLICING_BASE,"splicing","","base name for real subset",false);
    args.add_string(Opt::INTRONIC_BASE,"intronic","","base name for real subset",false);
    args.add_string(Opt::INTERGENIC_BASE,"intergenic","","base name for real subset",false);
    args.add_string(Opt::OUTPUT_FP,"output","","output file",true);

    if(argc <= 1 || strcmp(argv[1],"--help")==0){
        std::cerr<<args.get_help()<<std::endl;
        exit(1);
    }

    args.parse_args(argc,argv);

    // first create the execution string
    std::string cl="gtex_stats ";
    for (int i=0;i<argc;i++){
        if(i==0){
            cl+=argv[i];
        }
        else{
            cl+=" ";
            cl+=argv[i];
        }
    }

    QuantRes quantRes;
    if(args.is_set(Opt::REAL_BASE)){
        std::cout<<"loading real annotation"<<std::endl;
        std::string real_gff = args.get_string(Opt::REAL_BASE);
        real_gff = real_gff+".gtf";
        quantRes.load_gff(real_gff,false,true,false,false,false);
    }
    if(args.is_set(Opt::SPLICING_BASE)){
        std::cout<<"loading splicing annotation"<<std::endl;
        std::string splicing_gff = args.get_string(Opt::SPLICING_BASE);
        splicing_gff = splicing_gff+".gtf";
        quantRes.load_gff(splicing_gff,false,false,true,false,false);
    }
    if(args.is_set(Opt::INTRONIC_BASE)){
        std::cout<<"loading intronic annotation"<<std::endl;
        std::string intronic_gff = args.get_string(Opt::INTRONIC_BASE);
        intronic_gff = intronic_gff+".gtf";
        quantRes.load_gff(intronic_gff,false,false,false,true,false);
    }
    if(args.is_set(Opt::INTERGENIC_BASE)){
        std::cout<<"loading intergenic annotation"<<std::endl;
        std::string intergenic_gff = args.get_string(Opt::INTERGENIC_BASE);
        intergenic_gff = intergenic_gff+".gtf";
        quantRes.load_gff(intergenic_gff,false,false,false,false,true);
    }

    std::cout<<"loading base annotation"<<std::endl;
    quantRes.load_gff(args.get_string(Opt::ANNOTATION_FP),true,false,false,false,false);

    std::cout<<"loading reads"<<std::endl;
    quantRes.load_sim_reads(args.get_string(Opt::READS_FP));
    std::cout<<"loading stringtie assembly"<<std::endl;
    quantRes.load_strg(args.get_string(Opt::STRG_FP));
    std::cout<<"loading salmon abundances"<<std::endl;
    quantRes.load_slmn(args.get_string(Opt::SLMN_FP));
    std::cout<<"loading kallisto abundances"<<std::endl;
    quantRes.load_klst(args.get_string(Opt::KLST_FP));

    std::cout<<"computing TPMs"<<std::endl;
    quantRes.compute_tpm();
    std::cout<<"writing output"<<std::endl;
    quantRes.write(args.get_string(Opt::OUTPUT_FP));

    return 0;
}

// TODO: need to add another column which contains the number of reads from noise transcripts that overlap transcript in each row