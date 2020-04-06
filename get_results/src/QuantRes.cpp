//
// Created by sparrow on 3/8/20.
//

#include "QuantRes.h"
#include "FastaTools.h"

void QuantRes::get_poly_tid(std::string& recID,std::string& tid){
    // get read number which is first
    std::size_t prev_delim = 0;
    std::size_t delim = recID.find("/");

    // get transcript id
    prev_delim = delim;
    delim = recID.find(";", prev_delim + 1);
    tid = recID.substr(prev_delim + 1, (delim - prev_delim) - 1);
}

void QuantRes::load_sim_reads(std::string faFP){
    FastaReader fastaReader(faFP);
    bool first_read = true;
    while(fastaReader.good()){
        FastaRecord faRec;
        fastaReader.next(faRec);

        std::string tid;
        get_poly_tid(faRec.id_,tid);

        this->qr_it.first = this->quant_res.find(tid);
        if(this->qr_it.first == this->quant_res.end()){
            continue;
            std::cerr<<"unknown TID found: "<<tid<<std::endl;
            exit(-1);
        }
        else{
            if(first_read){
                this->readlen = faRec.seq_.size();
                first_read = false;
            }
            if(this->readlen != faRec.seq_.size()){
                std::cerr<<"reads of varying length not supported"<<std::endl;
                exit(-1);
            }
            this->qr_it.first->second.num_sim_reads++;
        }
    }
}

void QuantRes::load_gff(std::string gffFP,bool base, bool real, bool splicing, bool intronic, bool intergenic){
    FILE* gff_file = fopen(gffFP.c_str(), "r");
    if (gff_file == nullptr){
        std::cerr << "@ERROR::Couldn't open annotation: " << gffFP<< std::endl;
        exit(1);
    }
    GffReader gffReader;
    gffReader.init(gff_file,true);
    gffReader.readAll();

    GffObj *p_gffObj;

    for (int i = 0; i < gffReader.gflst.Count(); ++i) {
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count() == 0) {
            continue;
        }
        this->qr_it = this->quant_res.insert(std::make_pair(std::string(p_gffObj->getID()),QuantRec(p_gffObj, base, real, splicing, intronic, intergenic)));
        if(base){
            if(!this->qr_it.second){ // if already existed - mention that base exists
                this->qr_it.first->second.base = true;
            }
            this->qr_it.first->second.add_base(p_gffObj);
        }
        else{
            if(!this->qr_it.second){
                std::cerr<<"duplicate transcript IDs"<<this->qr_it.first->first<<std::endl;
                exit(-1);
            }
        }
    }
}

void QuantRes::load_slmn(std::string slmnFP){
    std::ifstream slmn_ss;
    slmn_ss.open(slmnFP.c_str(),std::ios::in);
    if (!slmn_ss.good()){
        std::cerr<<"@ERROR::Couldn't open file with the salmon quantification results: "<<slmnFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string line;

    std::getline(slmn_ss,line); // skip the header line

    while (std::getline(slmn_ss,line)) {
        // tid
        std::size_t found_tid = line.find("\t");
        if (found_tid==std::string::npos){std::cerr<<"ERROR #1 slmn"<<std::endl;exit(-1);}
        std::string tid = line.substr(0,found_tid);
        // length
        std::size_t found_len = line.find("\t",found_tid+1);
        if (found_len==std::string::npos){std::cerr<<"ERROR #2 slmn"<<std::endl;exit(-1);}
        // elen
        std::size_t found_elen = line.find("\t",found_len+1);
        if (found_elen==std::string::npos){std::cerr<<"ERROR #3 slmn"<<std::endl;exit(-1);}
        // tpm
        std::size_t found_tpm = line.find("\t",found_elen+1);
        if (found_tpm==std::string::npos){std::cerr<<"ERROR #4 slmn"<<std::endl;exit(-1);}
        std::string tpm_tmp = line.substr(found_elen+1,found_tpm - found_elen - 1);
        float tpm = std::atof(tpm_tmp.c_str());
        // num_reads
        std::string nr_tmp = line.substr(found_tpm+1,line.size() - 1);
        float nr = std::atof(nr_tmp.c_str());

        // find record
        this->qr_it.first = this->quant_res.find(tid);
        if(this->qr_it.first==this->quant_res.end()){
            std::cerr<<"salmon tid not found: "<<tid<<std::endl;
            exit(-1);
        }
        this->qr_it.first->second.slmn_reads = nr;
        this->qr_it.first->second.slmn_tpm = tpm;
    }
    slmn_ss.close();
}

void QuantRes::load_klst(std::string klstFP){
    std::ifstream klst_ss;
    klst_ss.open(klstFP.c_str(),std::ios::in);
    if (!klst_ss.good()){
        std::cerr<<"@ERROR::Couldn't open file with the kallisto quantification results: "<<klstFP<<std::endl;
        exit(1);
    }
    std::ios::sync_with_stdio(false);

    std::string line;

    std::getline(klst_ss,line); // skip the header line

    while (std::getline(klst_ss,line)) {
        // tid
        std::size_t found_tid = line.find("\t");
        if (found_tid==std::string::npos){std::cerr<<"ERROR #1 klst"<<std::endl;exit(-1);}
        std::string tid = line.substr(0,found_tid);
        // length
        std::size_t found_len = line.find("\t",found_tid+1);
        if (found_len==std::string::npos){std::cerr<<"ERROR #2 klst"<<std::endl;exit(-1);}
        // elen
        std::size_t found_elen = line.find("\t",found_len+1);
        if (found_elen==std::string::npos){std::cerr<<"ERROR #3 klst"<<std::endl;exit(-1);}
        // tpm
        std::size_t found_nr = line.find("\t",found_elen+1);
        if (found_nr==std::string::npos){std::cerr<<"ERROR #4 klst"<<std::endl;exit(-1);}
        std::string nr_tmp = line.substr(found_elen+1,found_nr - found_elen - 1);
        float nr = std::atof(nr_tmp.c_str());
        // num_reads
        std::string tpm_tmp = line.substr(found_nr+1,line.size() - 1);
        float tpm = std::atof(tpm_tmp.c_str());

        // find record
        this->qr_it.first = this->quant_res.find(tid);
        if(this->qr_it.first==this->quant_res.end()){
            std::cerr<<"kallisto tid not found: "<<tid<<std::endl;
            exit(-1);
        }
        this->qr_it.first->second.klst_reads = nr;
        this->qr_it.first->second.klst_tpm = tpm;
    }
    klst_ss.close();
}

void QuantRes::load_strg(std::string strgFP){
    FILE* strg_file = fopen(strgFP.c_str(), "r");
    if (strg_file == nullptr){
        std::cerr << "@ERROR::Couldn't open stringtie assembly results: " << strgFP<< std::endl;
        exit(1);
    }
    GffReader gffReader;
    gffReader.init(strg_file,true);
    gffReader.readAll(true);

    GffObj *p_gffObj;

    for (int i = 0; i < gffReader.gflst.Count(); ++i) {
        p_gffObj = gffReader.gflst.Get(i);
        if (p_gffObj->isDiscarded() || p_gffObj->exons.Count() == 0) {
            continue;
        }
        int ref_id_attid = p_gffObj->names->attrs.getId("reference_id");
        if(ref_id_attid==-1){ // atribute not found
            std::cerr<<"no attribute reference_id found in file"<<std::endl;
            exit(-1);
        }
        else{
            if(!p_gffObj->attrs->getAttr(ref_id_attid)){ // atribute not found - novel transcript
                continue;
            }

            std::string tid = p_gffObj->attrs->getAttr(ref_id_attid);

            this->qr_it.first = this->quant_res.find(tid);
            if(this->qr_it.first==this->quant_res.end()){
                std::cerr<<"reference transcript not found: "<<tid<<std::endl;
                exit(-1);
            }
            int len = p_gffObj->len();
            int elen = get_elen(p_gffObj);
            int cov_attid = p_gffObj->names->attrs.getId("cov");
            if(cov_attid==-1){std::cerr<<"coverage not found"<<std::endl;exit(-1);}
            float cov = std::atof(p_gffObj->attrs->getAttr(cov_attid));
            int tpm_attid = p_gffObj->names->attrs.getId("TPM");
            if(tpm_attid==-1){std::cerr<<"TPM not found"<<std::endl;exit(-1);}
            float tpm = std::atof(p_gffObj->attrs->getAttr(tpm_attid));
            this->qr_it.first->second.strg_len = len;
            this->qr_it.first->second.strg_elen = elen;
            this->qr_it.first->second.strg_cov = cov;
            this->qr_it.first->second.strg_tpm = tpm;
        }
    }
}

void QuantRes::write(std::string outFP){
    std::ofstream out_ss(outFP.c_str());
    out_ss<<"tid,len,elen_true,elen_base,sim_nreads,sim_tpm,strg_nreads,strg_tpm,slmn_nreads,slmn_tpm,klst_nreads,klst_tpm"<<std::endl;

    for(auto& rec : quant_res){
        if(rec.second.base || rec.second.real){
            out_ss << rec.second.get_rec_str(this->readlen)<<std::endl;
        }
    }

    out_ss.close();
}

void QuantRes::compute_tpm(){
    for(auto tx : this->quant_res){
        if(tx.second.real && !tx.second.base){
            std::cerr<<"real transcript found without a base: "<<tx.first<<std::endl;
            exit(-1);
        }
    }

    // first compute sum of RPKs
    float rpk_sum = 0;
    for(auto tx : this->quant_res){ // for each transcript
        if(tx.second.num_sim_reads>0){ // if reads were simulated for this transcripts whether real, intronic, splicing or intergenic
            rpk_sum = rpk_sum+tx.second.get_rpk();
        }
        else{ // otherwise we just skip
            continue;
        }
    }
    // next compute the per million scaling factor
    float pmf = rpk_sum/1000000.0;
    // next use pmf to get new TPM values
    for(auto tx : this->quant_res){ // for each transcript
        if(tx.second.num_sim_reads>0){ // if reads were simulated for this transcripts whether real, intronic, splicing or intergenic
            tx.second.set_sim_tpm(tx.second.get_rpk()/pmf);
        }
        else if(tx.second.base){ // if empty transcript but is part of the base annotation
            tx.second.set_sim_tpm(tx.second.get_rpk()/pmf);
        }
        else{ // otherwise we just skip
            continue;
        }
    }
}