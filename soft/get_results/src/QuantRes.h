//
// Created by Ales Varabyou on 3/8/20.
//

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <unordered_map>

#include "gff.h"
#include "GFaSeqGet.h"

#ifndef GET_RESULTS_QUANTRES_H
#define GET_RESULTS_QUANTRES_H

class QuantRec {
public:
    QuantRec(GffObj* pgo,bool base, bool real, bool splicing, bool intronic, bool intergenic){
        this->elen = get_elen(pgo);
        this->len = pgo->len();
        this->tid = pgo->getID();

        this->base = base;
        this-> real = real;
        this-> splicing = splicing;
        this-> intronic = intronic;
        this-> intergenic = intergenic;
    };
    ~QuantRec()=default;

    void set_sim_tpm(float tpm){
        this->sim_tpm = tpm;
    }

    float get_rpk(){
        return this->num_sim_reads/(this->elen/1000.0);
    }

    int base_elen = 0;
    void add_base(GffObj* pgo){
        this->base_elen = get_elen(pgo);
    }

    int num_sim_reads=0;
    float strg_cov=0;
    int strg_len=0;
    int strg_elen=0;
    float strg_tpm=0;
    float slmn_reads=0;
    float slmn_tpm=0;
    float klst_reads=0;
    float klst_tpm=0;

    float sim_tpm = 0; // simulated TPM

    bool base = false; // belongs to the base annotation
    bool real = false; // belongs to the real annotation
    bool splicing = false; // belongs to the splicing annotation
    bool intronic = false; // belongs to the intronic annotation
    bool intergenic = false; // belongs to the intergenic annotation

    std::string get_rec_str(int readlen){
        std::string rec_str = "";
        float num_strg_reads = ((float)strg_elen*strg_cov)/(float)readlen;
        if(strg_cov==0){
            num_strg_reads=0;
        }
        rec_str = tid+","+
                  std::to_string(len)+","+
                  std::to_string(elen)+","+
                  std::to_string(base_elen)+","+
                  std::to_string(num_sim_reads)+","+
                  std::to_string(sim_tpm)+","+
                  std::to_string(strg_elen)+","+
                  std::to_string(num_strg_reads)+","+
                  std::to_string(strg_tpm)+","+
                  std::to_string(slmn_reads)+","+
                  std::to_string(slmn_tpm)+","+
                  std::to_string(klst_reads)+","+
                  std::to_string(klst_tpm);
        return rec_str;
    }

private:
    std::string tid;
    int len = 0;
    int elen = 0;

    int get_elen(GffObj* pgo){
        int cur_len = 0;
        for(int i=0;i<pgo->exons.Count();i++){
            cur_len+=pgo->exons.Get(i)->len();
        }
        return cur_len;
    }
};

class QuantRes {
public:
    QuantRes()=default;
    ~QuantRes()=default;

    void load_sim_reads(std::string inFP);
    void load_gff(std::string gffFP,bool base, bool real, bool splicing, bool intronic, bool intergenic);
    void load_strg(std::string inFP);
    void load_slmn(std::string inFP);
    void load_klst(std::string inFP);

    void compute_tpm();

    void write(std::string outFP);

private:
    std::unordered_map<std::string,QuantRec> quant_res; // map of tid to quantification record
    std::pair<std::unordered_map<std::string,QuantRec>::iterator,bool> qr_it; // iterator to the above container

    int readlen = 0;

    void get_poly_tid(std::string& recID,std::string& tid);

    int get_elen(GffObj* pgo){
        int cur_len = 0;
        for(int i=0;i<pgo->exons.Count();i++){
            cur_len+=pgo->exons.Get(i)->len();
        }
        return cur_len;
    }
};


#endif //GET_RESULTS_QUANTRES_H
