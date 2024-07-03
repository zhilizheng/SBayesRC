/* SBayesRC
 *
 * This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
 * Develped by Zhili Zheng <zhilizheng@outlook.com>, 2021
 */

#include "SBayesRC.h"
#include "commR.h"

#include "dist.h"
#include <Eigen/Eigen>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <cmath>
#include "Timer.h"
#include <algorithm>


SBayesRC::SBayesRC(int niter, int burn, VectorXf fbhat, int numAnno, vector<string> &annoStrs, std::string mldmDir, double vary, VectorXf n, VectorXf fgamma, VectorXf pi, const std::vector<int> &rmSNPIndices, double starth2, double cutThresh, bool bOrigin, std::string outPrefix, std::string samVe, double resam_thresh, bool bOutDetail, int outFreq, double initAnnoSS){

    this->niter = niter;
    this->burn = burn;
    this->bOrigin = bOrigin;
    this->vary = vary;
    this->pi = pi;
    this->n = n;
    this->resam_thresh = resam_thresh;
    this->bOutDetail = bOutDetail;
    this->outPrefix = outPrefix;
    this->curSamVe = samVe;
    this->fgamma = fgamma;
    this->outFreq = outFreq;
    this->thinFreq = outFreq;
    this->b = fbhat;

    int nMarker = fbhat.size();
    ndist = fgamma.size();  // number of mixture distribution components

    //MatrixXf annoMat(nMarker, numAnno);
    bAnnot = false;
    if(numAnno > 0){
        bAnnot = true;
        Rcout << "Running SBayesRC with annotation" << std::endl;
    }else{
        Rcout << "Running SBayesRC without annotation" << std::endl;
    }

    string fileAnnot = outPrefix + ".annot.tmp.bin";

    if(bAnnot){
        snpPi.resize(nMarker, ndist);
        anno = new AnnoProb(fileAnnot, numAnno, pi, snpPi, bOutDetail, initAnnoSS);
        if(!outPrefix.empty()) anno->open(outPrefix, annoStrs);
    }
    //annoMat.resize(0, 0);

    // read LD
    blockLDeig.readLD(mldmDir, cutThresh, fbhat, "");
    m = blockLDeig.getNMarker();
    nBlocks = blockLDeig.getNBlocks();

    if(m != nMarker){
        Rcout << "inconsitent size between LD information and bhat" << std::endl;
        throw("Error");
    }

    varg = vary * starth2;
    Rcout << " Var_y: " << vary << endl;
    //float varg1 = vary * 0.005;
    //float cH2 = starth2 / nBlocks;
    //Rcout << " cH2: " << cH2 << endl;
    //vare.fill(vary * (1.0 - cH2));

    vare.resize(nBlocks);
    vare.fill(vary);
    Rcout << " Mean vare: " << vare.mean() << std::endl;


    sigmaSq = varg / (m * fgamma.dot(pi));
    //float sigmaSq = 0.01;
    //float sigmaSq = varg / (5045.0 * pi);

    //IntegerVector idx_dist = seq(0, ndist - 1);
    //Rcout << "idx_dist: " << idx_dist << std::endl;
    Rcout << "Number of distributions: " << ndist << std::endl;

    beta = VectorXf::Zero(m);
    //betasum = VectorXf::Zero(m);
    betasum_all = VectorXf::Zero(m);
    betamean = VectorXf::Zero(m);
    betaM2 = VectorXf::Zero(m);

    // history
    hsq_mcmc = VectorXf::Zero(niter);
    hsq2_mcmc = VectorXf::Zero(niter);
    pi_mcmc = MatrixXf::Zero(niter, ndist);
    n_comp_mcmc = MatrixXf::Zero(niter, ndist);
    vg_comp_mcmc = MatrixXf::Zero(niter, ndist);

    // iteration infos
    thinIter = (niter - burn) / thinFreq;
    iter_infos.resize(thinIter, 6);
    vare_infos.resize(thinIter, nBlocks);
    ssq_infos.resize(thinIter, nBlocks);
    hsq_infos.resize(thinIter, nBlocks);
    pip_count = MatrixXf::Zero(m, ndist);

    VectorXf abs_b = fbhat.array().abs();
    std::sort(abs_b.data(), abs_b.data() + abs_b.size());
    int index8 = 0.8 * (m - 1);
    betaThresh = abs_b[index8];
    Rcout << "Set beta cutoff threshold: " << betaThresh << std::endl;

    delSNPs.resize(nBlocks);
    if(rmSNPIndices.size() != 0){
        int curIdx = 0;
        int curRMIndex = rmSNPIndices[curIdx];
        bool b_allLoop = false;
        for(int blk = 0; blk < nBlocks; blk++){
            int start = blockLDeig.getIdxStart(blk);
            int end = blockLDeig.getIdxEnd(blk);
            if(b_allLoop) break;
            while(curRMIndex >= start && curRMIndex <= end){
                delSNPs[blk].insert(curRMIndex);
                curIdx++;
                if(curIdx < rmSNPIndices.size()){
                    curRMIndex = rmSNPIndices[curIdx];
                }else{
                    b_allLoop = true;
                    break;
                }
            }
        }

        Rcout << curIdx << " SNPs will be removed from MCMC." << std::endl;
    }
}

const int index_size = sizeof(uint32_t);
const int value_size = sizeof(float);

void SBayesRC::setOutBeta(bool bOut){
    this->bOutBeta = bOut;
    if(bOutBeta){
        this->fpBeta = fopen((outPrefix + ".beta.bin").c_str(), "wb");
        if(!fpBeta){
            Rcout << "Error: can't open [" << outPrefix << ".beta.bin] for writing" << std::endl;
            throw("ERROR");
        }
        uint32_t size[2] = {thinIter, m};
        fwrite(size, sizeof(uint32_t), 2, fpBeta);

        head_byte = 2 * index_size;
        content_byte = head_byte + thinIter * index_size;
        fseek(fpBeta, content_byte, SEEK_SET);
    }
}

void SBayesRC::writeBeta(const VectorXf& beta, int curIter){
    int curIterThin = (curIter - burn) / thinFreq;

    Eigen::SparseVector<float> sparse_vec = beta.sparseView();
    uint32_t curRow = sparse_vec.nonZeros();
    if(fwrite(sparse_vec.innerIndexPtr(), index_size, curRow, fpBeta) != curRow){
        Rcout << "Error: can't write index to MCMC beta values" << std::endl;
        throw("ERROR");
    }
    if(fwrite(sparse_vec.valuePtr(), value_size, curRow, fpBeta) != curRow){
        Rcout << "Error: can't write value to MCMC beta values" << std::endl;
        throw("ERROR");
    }
    long curPos = ftell(fpBeta);

    fseek(fpBeta, head_byte + curIterThin * index_size, SEEK_SET);
    fwrite(&curRow, sizeof(uint32_t), 1, fpBeta);

    fseek(fpBeta, curPos, SEEK_SET);

}

void SBayesRC::mcmc(){
    Rcout << "Start MCMC with " << niter << " iterations..." << endl;

    bool estimatePi = true;

    float nub = 4;
    float nue = 4;
    float scaleb = (nub - 2) / nub * sigmaSq;
    VectorXf scalee = (nue - 2) / nue * vare.array();

    Timer timer;
    timer.start("sbrc");

    vector<VectorXf> whats;
    for(int i = 0; i < nBlocks; i++){
        VectorXf temp = VectorXf::Zero(blockLDeig.getq(i));
        whats.emplace_back(temp);
    }

    MatrixXf cur_causal = MatrixXf::Zero(m, ndist);
    MatrixXf z = MatrixXf::Zero(m, ndist - 1);
    MatrixXf vg_snp_comp = MatrixXf::Zero(m, ndist);
    VectorXf vg_comp(ndist); 
    VectorXf n_comp(ndist);
    VectorXf Ve_param(nBlocks);

    
    for(int iter = 0; iter < niter; iter++){
        //Rcout << "------------------------------------" << std::endl;
        //Rcout << "iterStart " << iter << std::endl;
        //MatrixXf z_dist = MatrixXf::Zero(m, ndist - 1);
        cur_causal.setZero();
        z.setZero();
        //int nnz = 0;

        //VectorXf numSnpDist = VectorXf::Zero(ndist);
        //VectorXf logPi = pi.array().log();
        MatrixXf logPiA;
        VectorXf logPi;
        if(bAnnot){
            logPiA = snpPi.array().log();
            //logPi = pi.array().log();
        }else{
            logPi = pi.array().log();
        }

        //VectorXf logPiComp = (1.0 - pi.array()).log(); 
        VectorXf wtSigmaSq;
        if(bOrigin){
            wtSigmaSq = fgamma.array() * 0.01 * varg;
        }else{
            wtSigmaSq = fgamma.array() * sigmaSq;
        }
        VectorXf invSigmaSq = wtSigmaSq.array().inverse();
        VectorXf logSigmaSq = wtSigmaSq.array().log();

        //Rcout << "  iter: " << iter <<  ", DEBUG clean memory" << std::endl; 
        #pragma omp parallel for schedule(dynamic)
        for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
            whats[idxBlk].setZero();
            //single_whats[idxBlk].setZero();
        }

        //Rcout << "  prepare: " << (timer.now() - tic)/1e9 << std::endl;

        //Rcout << "  iter: " << iter <<  ", DEBUG start main" << std::endl; 
        #pragma omp parallel for schedule(dynamic)
        for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
            //Rcout << "Blk: " << idxBlk << std::endl;
            Ref<const MatrixXf> curQ = blockLDeig.getQ(idxBlk);
            Ref<VectorXf> wcorr = blockLDeig.getW(idxBlk);
            std::set<int>& delSNP = delSNPs.at(idxBlk);
            int startM = blockLDeig.getIdxStart(idxBlk);
            int endM = blockLDeig.getIdxEnd(idxBlk);

            int baseStart = blockLDeig.getStartPos(idxBlk);
            float q = blockLDeig.getq(idxBlk);
            //int curNnz = 0;

            Ref<VectorXf> what = whats[idxBlk];
            float vareDn = n[idxBlk] / vare[idxBlk];
            VectorXf invLhs = 1.0/(vareDn + invSigmaSq.array());
            VectorXf logInvLhsMsigma = invLhs.array().log() - logSigmaSq.array();
            //VectorXf sqrtInvLhs = invLhs.array().sqrt();

            for(int idx = startM; idx <= endM; idx++){
                int found = delSNP.count(idx);
                if(found != 0){
                    continue;
                }
                float oldSample = beta[idx];
                Ref<const VectorXf> curQi = curQ.col(idx - baseStart);
                float rhs = (curQi.dot(wcorr) + oldSample)*vareDn;
                VectorXf uhat = invLhs.array() * rhs;
                VectorXf logDelta;
                if(bAnnot){
                    logDelta = 0.5*(logInvLhsMsigma.array() + uhat.array() * rhs) + logPiA.row(idx).transpose().array();
                    //Rcout << "iter: " << iter << ", idx " << idx << ": " << logPiA.row(idx) << std::endl;
                    //Rcout << "    " << logPi.transpose() << std::endl;
                    logDelta[0] = logPiA(idx, 0);
                    //VectorXf logDelta2 = 0.5*(logInvLhsMsigma.array() + uhat.array() * rhs) + logPi.array();
                    //logDelta2[0] = logPi(0);
                    //Rcout << "anno: " << logDelta.transpose() << std::endl;
                    //Rcout << "noan: " << logDelta2.transpose() << std::endl;
                }else{
                    logDelta = 0.5*(logInvLhsMsigma.array() + uhat.array() * rhs) + logPi.array();
                    logDelta[0] = logPi(0);
                }

                //VectorXf logDelta = 0.5*(logInvLhsMsigma.array() + uhat.array() * rhs) + logPi.array();
                //logDelta[0] = logPi(idx, 0);
                //NumericVector probDelta(ndist);
                VectorXf probDelta(ndist);

                for(int k = 0; k < ndist; k++){
                    probDelta[k] = 1.0/(logDelta.array() - logDelta[k]).exp().sum();
                    /*
                    if(isnan(probDelta[k])){
                        Rcout << " iter: " << iter << ", block: " << idxBlk << ", idx: " << idx << ", dist: " << k << ", nan" << std::endl;
                    }
                    */
                }
                int delta = Bernoulli::sample(probDelta);

                //numSnpDist[delta]++;
                cur_causal(idx, delta) = 1;

                //stop("check output");
                float adj_wcorr = 0;
                bool bReviseWhat = false;
                if(delta != 0){
                    beta[idx] = Normal::sample(uhat[delta], invLhs[delta]);

                    /*
                    float rate_b = abs((beta[idx] - b[idx])/b[idx]);
                    bool sameSign = (b[idx] >= 0.0f) == (beta[idx] >= 0.0f);
                    float compare_rate = sameSign ? 2 : 1.1;

                    if(abs(beta[idx]) > betaThresh && rate_b > compare_rate){
                        adj_wcorr = oldSample;
                        beta[idx] = 0;
                        cur_causal(idx, delta) = 0;
                        delSNP.insert(idx);
                    }else{
                    */
                        //wcorr = wcorr + curQi * (oldSample - beta[idx]);
                        adj_wcorr = oldSample - beta[idx];
                        //vg_snps(idx) = beta[idx] * beta[idx] * weightQ(idx);
                        //vg_snps2(idx) = beta[idx] * beta[idx] * weightQ2(idx);
                        what = what + curQi * beta[idx];
                        bReviseWhat = true;
                        //curNnz = curNnz + 1;

                        //z(idx, 0) = 1;
                        //if(delta > 1) z(idx, 1) = 1;
                        //if(delta > 2) z(idx, 2) = 1;
                        for(int i2 = 0; i2 < delta ; i2++){
                            z(idx, i2) = 1;
                        }
                    //}
                }else{
                    //if(oldSample != 0) wcorr = wcorr + curQi * oldSample;
                    adj_wcorr = oldSample;
                    beta[idx] = 0;
                    //vg_snps(idx) = 0;
                    //vg_snps2(idx) = 0;
                }
                if(adj_wcorr != 0){
                    wcorr = wcorr + curQi * adj_wcorr;
                }
           }

       }

        VectorXf Vg_block(nBlocks);
        #pragma omp parallel for schedule(dynamic)
        for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
            float curVarg = whats[idxBlk].dot(whats[idxBlk]);
            Vg_block[idxBlk] = curVarg;
        }

        //if(iter >= 20){
            betasum_all = betasum_all + beta;
        //}
        if(iter >= burn) {
            //betasum = betasum + beta;
            VectorXf delta = beta - betamean;
            betamean = betamean.array() + delta.array() / (iter - burn + 1);
            betaM2 = betaM2.array() + delta.array() * (beta.array() - betamean.array());
            pip_count = pip_count + cur_causal;
            if(bOutBeta){
                if(!((iter+1) % thinFreq)){
                    writeBeta(beta, iter);
                }
            }
        }

        //if(estimateSigmaSq) sigmaSq = (dot(beta, beta) + nub*scaleb)/rchisq(1, nnz+nub)[0];
        //if(estimateSigmaSq) sigmaSq = median(cur_SigmaSqs);
        // float hsq = 0.5 * cur_SigmaSqs.sum();
        //Rcout << "  iter: " << iter <<  ", DEBUG V" << std::endl; 
        varg = Vg_block.sum();
        VectorXf vg_snps = beta.array().square();
        float varg2 = vg_snps.sum();
        float hsq2 = varg2 / vary;

        // ssq
        VectorXf ssq_block(nBlocks);
        #pragma omp parallel for schedule(dynamic)
        for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
            int curStart = blockLDeig.getIdxStart(idxBlk);
            int curEnd = blockLDeig.getIdxEnd(idxBlk);
            ssq_block[idxBlk] = vg_snps.segment(curStart, curEnd - curStart + 1).sum();
        }
 
        // resample Ve
        if(curSamVe == "fixVe"){
            for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                Ve_param[idxBlk] = -1;
            }
        }else if(curSamVe == "samVe"){
            #pragma omp parallel for schedule(dynamic)
            for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                Ref<VectorXf> wcorr = blockLDeig.getW(idxBlk);
                Ve_param[idxBlk] = n[idxBlk] * wcorr.dot(wcorr) + nue * scalee[idxBlk];
            }
        }else if(curSamVe == "allMixVe"){
            #pragma omp parallel for schedule(dynamic)
            for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                if(Vg_block[idxBlk] < 1e-8 || ssq_block[idxBlk] < 1e-8){
                    Ve_param[idxBlk] = -1;
                }else if(ssq_block[idxBlk] / Vg_block[idxBlk] > resam_thresh){
                    Ref<VectorXf> wcorr = blockLDeig.getW(idxBlk);
                    Ve_param[idxBlk] = n[idxBlk] * wcorr.dot(wcorr) + nue * scalee[idxBlk];
                }else{
                    Ve_param[idxBlk] = -1;
                }
            }
        }else if(curSamVe == "mixVe"){
            const int tryIter = 50;
            if(iter < tryIter){
                for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                    Ve_param[idxBlk] = -1;
                }
            }else if(iter == tryIter){
                int nSam = 0;
                for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                    if(Vg_block[idxBlk] < 1e-8 || ssq_block[idxBlk] < 1e-8){
                        Ve_param[idxBlk] = -1;
                    }else if(ssq_block[idxBlk] / Vg_block[idxBlk] > resam_thresh){
                        Ve_param[idxBlk] = 1;
                        nSam++;
                    }else{
                        Ve_param[idxBlk] = -1;
                    }
                }

                Rcout << "After " << tryIter << " iterations' trying, " << nSam << " blocks need to be controlled by sampling of Ve. " << std::endl;
            }

            if(iter >= tryIter){
                #pragma omp parallel for schedule(dynamic)
                for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                    if(Ve_param[idxBlk] > 0){
                        Ref<VectorXf> wcorr = blockLDeig.getW(idxBlk);
                        Ve_param[idxBlk] = n[idxBlk] * wcorr.dot(wcorr) + nue * scalee[idxBlk];
                    }
                }
            }
        }else{
            Rcout << "Error sample Ve manner" << std::endl;
            throw("Error");
        }

        for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
            float curParam = Ve_param[idxBlk];
            if(curParam > 0){
                float sample1 = InvChiSq::sample(blockLDeig.getq(idxBlk) + nue, curParam);
                if(sample1 / vary > 0.7){
                    vare[idxBlk] = sample1; 
                }else{
                    vare[idxBlk] = vary;
                }
            }else{
                vare[idxBlk] = vary;
            }
        }


        /*
        MatrixXf vg_snp_comp = cur_causal.array().colwise() * vg_snps.array();
        VectorXf vg_comp = vg_snp_comp.colwise().sum().array() / varg2; 
        VectorXf n_comp = cur_causal.colwise().sum();
        */
        #pragma omp parallel for
        for(int i = 0; i < ndist; i++){
            vg_snp_comp.col(i) = cur_causal.col(i).array() * vg_snps.array();
            vg_comp(i) = vg_snp_comp.col(i).sum() / varg2;
            n_comp(i) = cur_causal.col(i).sum();
        }

        //Rcout << "  finish Vg cal: " << (timer.now() - tic)/1e9 << std::endl;
        //Rcout << "  iter: " << iter <<  ", DEBUG sigmasq" << std::endl; 
        float hsq = varg / vary;
        //estimateSigmaSq = false;
        if(estimateSigmaSq){
            //sigmaSq = varg / (m * fgamma.dot(pi));
            //sigmaSq = varg / (snpPi.array().rowwise() * pi.transpose().array()).sum();
            /*
            if(bAnnot){
                sigmaSq = varg / fgamma.dot(snpPi.colwise().sum());
            }else{
                sigmaSq = varg / (m * fgamma.dot(pi));
            }
            */
            sigmaSq = varg / (m * fgamma.dot(pi));
        }

        //Rcout << "  iter: " << iter <<  ", DEBUG Pi" << std::endl; 
        //estimatePi = false;
        if(estimatePi){
            if(bAnnot){
                anno->samplePi(z, snpPi);
                #pragma omp parallel for
                for(int i = 0; i < ndist; i++){
                    pi(i) = snpPi.col(i).sum() /  m;
                }
                //pi = snpPi.colwise().sum().array() / m;
            }else{
                VectorXf n_comp_1 = n_comp.array() + 1;
                pi = Dirichlet::sample(ndist, n_comp_1);
                /*
                VectorXf ps(ndist);
                float sx = 0.0;
                for (int i = 0; i < ndist; i++){
                    ps[i] = rgamma(1, numSnpDist[i]+1, 1.0)[0];
                    sx += ps[i];
                }
                pi = ps / sx;
                */
            }
       }
       //Rcout << "  finish Pi: " << (timer.now() - tic)/1e9 << std::endl;

        hsq_mcmc[iter] = hsq;
        hsq2_mcmc[iter] = hsq2;
        pi_mcmc.row(iter) = pi;
        n_comp_mcmc.row(iter) = n_comp;
        vg_comp_mcmc.row(iter) = vg_comp;

        //Rcout << " Total time " << (timer.now() - tic)/1e9 << std::endl;
        //tic = timer.now();
        if(!((iter+1) % outFreq)){
            float m_vare = vare.mean();
            int nnz = n_comp.sum() - n_comp[0];
            if(!((iter+1) % (outFreq * 10))){
                double t100 = timer.elapse("sbrc");
                //Rprintf("\n iter %i, pi = %6.3f, nnz = %i, sigmaSq = %6.3f, hsq = %6.3f, vare = %6.3f, varg = %6.3f, time = %6.3f\n", iter + 1, pi, nnz, sigmaSq, hsq, m_vare, varg, t100); 
                //Rprintf("\n iter %i, nnz = %i, sigmaSq = %6.3f, hsq = %6.3f, vare = %6.3f, varg = %6.3f,  varg2 = %6.3f, time = %6.3f\n", iter + 1, nnz, sigmaSq, hsq, m_vare, varg, varg2, t100); 
                // check the SNPs
                int NdelSNPs = 0;
                VectorXf betaVals = betasum_all.array() / (iter + 1);
                VectorXi vDelSNPs = VectorXi::Zero(nBlocks);
                VectorXi signDelSNPs = VectorXi::Zero(nBlocks);
                float rate_thresh1 = 0, rate_thresh2 = 0;

                if(iter < 300){
                    rate_thresh1 = 4.0;
                    rate_thresh2 = 2.0;
                }else if(iter < 600){
                    rate_thresh1 = 3.0;
                    rate_thresh2 = 1.5;
                }else if(iter < 900){
                    rate_thresh1 = 2.0;
                    rate_thresh2 = 1.3;
                }else{
                    rate_thresh1 = 1.5;
                    rate_thresh2 = 1.1;
                }
                //rate_thresh1 = 4.0;
                //rate_thresh2 = 1.5;

                #pragma omp parallel for schedule(dynamic)
                for(int idxBlk = 0; idxBlk < nBlocks; idxBlk++){
                    Ref<const MatrixXf> curQ = blockLDeig.getQ(idxBlk);
                    Ref<VectorXf> wcorr = blockLDeig.getW(idxBlk);
                    std::set<int>& delSNP = delSNPs.at(idxBlk);
                    int startM = blockLDeig.getIdxStart(idxBlk);
                    int endM = blockLDeig.getIdxEnd(idxBlk);

                    int baseStart = blockLDeig.getStartPos(idxBlk);

                    for(int idx = startM; idx <= endM; idx++){
                        int found = delSNP.count(idx);
                        if(found != 0){
                            vDelSNPs[idxBlk] += 1;
                            continue;
                        }
                        if((iter + 1) == niter){
                            continue;
                        }
                        Ref<const VectorXf> curQi = curQ.col(idx - baseStart);
                        float betaVal = betaVals[idx];
                        float rate_b = abs((betaVal - b[idx])/b[idx]);
                        bool sameSign = (b[idx] >= 0.0f) == (betaVal >= 0.0f);
                        float compare_rate = sameSign ? rate_thresh1 : rate_thresh2;

                        if(abs(betaVal) > betaThresh && rate_b > compare_rate){
                            //Rcout << "DEL:" << "\t" << betaVal << "\t" << b[idx] << "\t" << rate_b << "\t" << compare_rate << std::endl;
                            float adj_wcorr = beta[idx];
                            wcorr = wcorr + curQi * adj_wcorr;
                            beta[idx] = 0;
                            betamean[idx] = 0;
                            betaM2[idx] = 0;
                            pip_count.row(idx).setZero();
                            betasum_all[idx] = 0;
                            delSNP.insert(idx);
                            vDelSNPs[idxBlk] += 1;
                            if(!sameSign) signDelSNPs[idxBlk] += 1;
                        }
                    }
                }

 
                //
                string n_str = "";
                string vg_str = "";
                for(int i = 1; i < ndist; i++){
                    n_str = n_str + "n" + to_string(i+1) + "=" + to_string((int)n_comp[i]) + ", ";
                    string vg_temp0 = to_string(vg_comp[i]);
                    string vg_temp1 = vg_temp0.substr(0, vg_temp0.find(".") + 3 + 1);
                    vg_str = vg_str + "vg" + to_string(i+1) + "=" + vg_temp1 + ", ";
                }
                Rprintf("  Iter %i, nnz=%i, sigmaSq=%.3f, hsq=%.3f, ssq=%.3f, %s%s vare=%.3f, rmVariants=%i, time=%.3f\n", iter + 1, nnz, sigmaSq, hsq, hsq2, n_str.c_str(), vg_str.c_str(), m_vare, vDelSNPs.sum(), t100); 

                //Rcout << "ProbDelta Mean: " << probDeltas.mean() << std::endl;
                timer.start("sbrc");
            }
            if(iter >= burn){
                int curRow = (iter - burn) / outFreq;
                iter_infos(curRow, 0) = hsq;
                iter_infos(curRow, 1) = nnz;
                iter_infos(curRow, 2) = sigmaSq;
                iter_infos(curRow, 3) = hsq2;
                iter_infos(curRow, 4) = m_vare;
                iter_infos(curRow, 5) = varg;

                vare_infos.row(curRow) = vare.cast<double>();
                hsq_infos.row(curRow) = Vg_block.cast<double>();
                ssq_infos.row(curRow) = ssq_block.cast<double>();

                //output anno
                if(bAnnot && (!outPrefix.empty())){
                    anno->computeProb();
                    if(bOutDetail) {
                        anno->computeDist(z, n_comp);
                    }
                    anno->computeEnrichBin(vg_snp_comp, varg2);
                    anno->computeEnrichQt(vg_snps, varg2);
                    anno->output();
                }
           }
        }
   }
    //vec mean_par = conv_to<vec>::from((sum(iter_infos, 0) / iter_infos.n_rows).as_col());
    //vec betaMean = conv_to<vec>::from(sum(beta_mcmc, 1) / beta_mcmc.n_cols);
    //NumericVector mean_par = as<NumericVector>(wrap((sum(iter_infos, 0) / iter_infos.n_rows).as_col()));
   // clean
   if(anno){
       annoSS = anno->getSigmaSq();
       if(!outPrefix.empty())anno->close();
       delete anno;
   }

   if(bOutBeta){
       fclose(fpBeta);
   }
}

VectorXd SBayesRC::get_mean_par_vec(){
    return  iter_infos.colwise().mean();
}

VectorXd SBayesRC::get_sd_par_vec(){
    VectorXd colMean = get_mean_par_vec();
    MatrixXd centered = iter_infos.rowwise() - colMean.transpose();
    VectorXd colVariance = (centered.array().square().colwise().sum()) / (iter_infos.rows() - 1);
    return colVariance.array().sqrt();
}


VectorXf SBayesRC::get_betaMean_vec(){
    return betasum.array() / (niter - burn);
}

VectorXf SBayesRC::get_betaMean_dir(){
    return betamean;
}

VectorXf SBayesRC::get_betaSD_dir(){
    return (betaM2.array() / (niter - burn - 1)).sqrt();
}

VectorXf SBayesRC::get_betaMean2_vec(){
    return betasum2.array() / (thinIter);
}


VectorXf SBayesRC::get_betaMean3_vec(){
    return betasum3.array() / (thinIter * 2);
}

VectorXf SBayesRC::get_hsq_mcmc(){
    return hsq_mcmc;
}

VectorXf SBayesRC::get_hsq2_mcmc(){
    return hsq2_mcmc;
}

VectorXf SBayesRC::get_beta(){
    return beta;
}

MatrixXf SBayesRC::get_pi_mcmc(){
    return pi_mcmc;
}

MatrixXf SBayesRC::get_n_comp_mcmc(){
    return n_comp_mcmc;
}

MatrixXf SBayesRC::get_vg_comp_mcmc(){
    return vg_comp_mcmc;
}

MatrixXf SBayesRC::get_pip(){
    return pip_count.array() / (niter - burn);
}

MatrixXf SBayesRC::get_pip2(){
    return pip_count2.array() / (thinIter);
}

MatrixXf SBayesRC::get_pip3(){
    return pip_count3.array() / (thinIter * 2);
}

MatrixXd SBayesRC::get_vare_infos(){
    return vare_infos;
}

MatrixXd SBayesRC::get_ssq_infos(){
    return ssq_infos;
}

MatrixXd SBayesRC::get_hsq_infos(){
    return hsq_infos;
}

VectorXd SBayesRC::get_anno_ss(){
    return annoSS;
}

void SBayesRC::outRmIndex(std::string outPrefix){
    ofstream outs((outPrefix + ".rm.snpidx").c_str());
    for(int idx = 0; idx < nBlocks; idx++){
        for(int index : delSNPs[idx]){
            outs << index << std::endl;
        }
    }
    outs.close();
}
