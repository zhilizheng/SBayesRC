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

#include "commR.h"
#include "dist.h"
#include "AnnoProb.h"
#include <cmath>
#include <boost/algorithm/string/join.hpp>

void AnnoProb::initP_Pi_annoAlpha(const VectorXf &pi, MatrixXf &snpPi) {
    int ndist = pi.size();
    float pi_sum = pi.sum();
    if (std::abs(pi_sum - 1) > 1e-4) {
        throw("Error: the sum of pi is not 1");
    }

    for (int i = 0; i < ndist; i++) {
        snpPi.col(i).setConstant(pi[i]);
    }

    VectorXf p(numComp);
    for (int i = 0; i < numComp; i++) {
        p[i] = pi.tail(ndist - 1 - i).sum() / pi.tail(ndist - i).sum();
    }

    for (int i = 0; i < numComp; i++) {
        snpP.col(i).setConstant(p[i]);
        alpha(0, i) = Normal::quantile_01(p[i]);
    }
}

void AnnoProb::annoCondProb_probit() {
    for (int i = 0; i < numComp; i++) {
        const VectorXf &alphai = alpha.col(i);
        annoCondProb(0, i) = Normal::cdf_01(alphai[0]);
        #pragma omp parallel for
        for (int j = 1; j < numAnno; j++) {
            annoCondProb(j, i) = Normal::cdf_01(alphai[0] + annoSD[j] * alphai[j]);
        }
    }
}

void AnnoProb::output(){
    // output alpha
    // todo: output annoCondProb, annoJointProb
    outs["vg_tot_annot"] << vg_tot_annot.transpose() << std::endl;
    //outs["vg_enrich_bin"] << vg_enrich_bin.transpose() << std::endl;
    outs["vg_enrich_qt"] << vg_enrich_qt.transpose() << std::endl;

    
    for(int i = 0; i < numDist; i++){
        if(bOutDetail){
            outs["annoDists" + to_string(i)] << annoDists.col(i).transpose() <<  std::endl;
        }
        outs["vg_annot_comp" + to_string(i)] << vg_annot_comp.col(i).transpose() << std::endl;
        outs["annoJointProb" + to_string(i)] << annoJointProb.col(i).transpose() << std::endl;
    }

    for(int i = 0; i < numComp; i++){
        outs["alpha" + to_string(i)] << alpha.col(i).transpose() << std::endl;
        outs["annoCondProb" + to_string(i)] << annoCondProb.col(i).transpose() << std::endl;
    }

    
}

void AnnoProb::open(string prefix, const vector<string> &annoStrs){
    // vg_tot_annot: numAnno * 1
    // vg_enrich_qt: numAnno * 1 
    // vg_enrich_bin: numAnno * 1
    //
    // annoDists: numAnno * numDist
    // vg_annot_comp: numAnno * numDist
    //
    // alpha: numAnno * numComp
    // annoCondProb: numAnno * numComp
    // annoJointProb: numAnno * numDist
    //
    outs.emplace("vg_tot_annot", ofstream((prefix + ".mcmcsamples.AnnoTotalGenVar").c_str()));
    //outs.emplace("vg_enrich_bin", ofstream((prefix + ".vg.enrich.bin").c_str()));
    outs.emplace("vg_enrich_qt", ofstream((prefix + ".mcmcsamples.AnnoPerSnpHsqEnrichment").c_str()));

    for(int i = 0; i < numDist; i++){
        if(bOutDetail){
            outs.emplace("annoDists" + to_string(i), ofstream((prefix + ".annoDists" + to_string(i+1) )));
        }
        outs.emplace("vg_annot_comp" + to_string(i), ofstream((prefix + ".mcmcsamples.AnnoGenVar_pi" + to_string(i+1) )));
        outs.emplace("annoJointProb" + to_string(i), ofstream((prefix + ".mcmcsamples.AnnoJointProb_pi" + to_string(i+1) )));
    }
 
    for(int i = 0; i < numComp; i++){
        outs.emplace("alpha" + to_string(i), ofstream((prefix + ".mcmcsamples.AnnoEffects_p" + to_string(i+1) )));
        outs.emplace("annoCondProb" + to_string(i), ofstream((prefix + ".mcmcsamples.AnnoCondProb_p" + to_string(i+1) )));
    }

    writeHeader(annoStrs);
}

void AnnoProb::writeHeader(const vector<string> &annoStrs){
    string delim = " ";
    string header = boost::algorithm::join(annoStrs, delim);

    outs["vg_tot_annot"] << header << std::endl;
    //outs["vg_enrich_bin"] << header << std::endl;
    outs["vg_enrich_qt"] << header << std::endl;

    for(int i = 0; i < numDist; i++){
        if(bOutDetail){
            outs["annoDists" + to_string(i)] << header << std::endl;
        }
        outs["vg_annot_comp" + to_string(i)] << header << std::endl;
        outs["annoJointProb" + to_string(i)] << header << std::endl;
    }

    for(int i = 0; i < numComp; i++){
        outs["alpha" + to_string(i)] << header << std::endl;
        outs["annoCondProb" + to_string(i)] << header << std::endl;
    }
}

void AnnoProb::close(){
    outs["vg_tot_annot"].close();
    //outs["vg_enrich_bin"].close();
    outs["vg_enrich_qt"].close();

    for(int i = 0; i < numDist; i++){
        if(bOutDetail){
            outs["annoDists" + to_string(i)].close();
        }
        outs["vg_annot_comp" + to_string(i)].close();
        outs["annoJointProb" + to_string(i)].close();
    }

    for(int i = 0; i < numComp; i++){
        outs["alpha" + to_string(i)].close();
        outs["annoCondProb" + to_string(i)].close();
    }
}

void AnnoProb::annoSigmaSS_sample() {
    // init sigmaSS: SetOnes(ssq.size());
    const int df = 4;
    const int scale = 1; // changed from 1

    float dfTilde = df + (numAnno - 1);  // exclude the intercept
    VectorXf scaleTilde = ssq.array() + df * scale;
    // #pragma omp parallel for
    for (int i = 0; i < numComp; i++) {
        //float scaleTilde = ssq[i] + df*scale;
        //values[i] = InvChiSq::sample(dfTilde, scaleTilde);
        //sigmaSq[i] = scaleTilde[i] / rchisq(1, dfTilde)[0];
        sigmaSq[i] = InvChiSq::sample(dfTilde, scaleTilde[i]);
    }
}

void AnnoProb::annoEffect_sample_Gibbs(const MatrixXf &z) {
    // alpha: numAnno * numComp
    // ssq: numComp
    // snpP: numSnps * numComp

    VectorXf numOnes_temp = z.colwise().sum();
    VectorXf numOnes(numComp + 1);
    numOnes << numSnps, numOnes_temp;

    // #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numComp; i++) {
        //VectorXf &alphai = alpha.col(i);
        Ref<VectorXf> alphai = alpha.col(i);
        //MatrixXf annoMati;
        MatrixXf *annotMatP;
        MatrixXf annoMatPO;

        int numDP = numOnes[i];

        //annoMati.setZero(numDP, numAnno);
        //zi.setZero(numDP);

        // get annotation coefficient matrix for component i
        if(numDP != 0){
            VectorXf y, zi;
            y.setZero(numDP);
            if (i == 0) {
                //annoMati = annoMat;
                annotMatP = &annoMat;
                zi = z.col(i);
            } else {
                annoMatPO.setZero(numDP, numAnno);
                zi.setZero(numDP);
                for (int j = 0, idx = 0; j < numSnps; ++j) {
                    if (z(j, i - 1)) {
                        annoMatPO.row(idx) = annoMat.row(j);
                        zi[idx] = z(j, i);
                        ++idx;
                    }
                }
                annotMatP = &annoMatPO;
            }
            MatrixXf &annoMati = (*annotMatP);

            VectorXf annoDiagi(numAnno);
            annoDiagi[0] = numDP;
            #pragma omp parallel for
            for (int k = 1; k < numAnno; ++k) {
                annoDiagi[k] = annoMati.col(k).squaredNorm();
            }

            // compute the mean of truncated normal distribution
            VectorXf mean = annoMati * alphai;

            // sample latent variables
            // # pragma omp parallel for schedule(dynamic)
            #pragma omp parallel for
            for (int j = 0; j < numDP; ++j) {
                if (zi[j])
                    y[j] = TruncatedNormal::sample_lower_truncated(mean[j], 1.0, 0.0);
                else
                    y[j] = TruncatedNormal::sample_upper_truncated(mean[j], 1.0, 0.0);
            }

            // adjust the latent variable by all annotation effects;
            y -= mean;

            // intercept is fitted with a flat prior
            float oldSample = alphai[0];
            float rhs = y.sum() + annoDiagi[0] * oldSample;
            float invLhs = 1.0 / annoDiagi[0];
            float ahat = invLhs * rhs;
            alphai[0] = Normal::sample(ahat, invLhs);
            y.array() += oldSample - alphai[0];
            //        cout << i << " alphai[0] " << alphai[0] << endl;

            // annotations are fitted with a normal prior
            ssq[i] = 0;
            for (int k = 1; k < numAnno; ++k) {
                oldSample = alphai[k];
                rhs = annoMati.col(k).dot(y) + annoDiagi[k] * oldSample;
                invLhs = 1.0 / (annoDiagi[k] + 1.0 / sigmaSq[i]);
                ahat = invLhs * rhs;
                alphai[k] = Normal::sample(ahat, invLhs);
                y += annoMati.col(k) * (oldSample - alphai[k]);
                ssq[i] += alphai[k] * alphai[k];
            }
        }else{
            alphai.setZero();
            alphai[0] = -10.0;
            ssq[i] = 0;
        }

        //cout << i << " " << alphai.transpose() << endl;

        /*
        for (int j = 0; j < numSnps; ++j) {
            snpP(j, i) = Normal::cdf_01(annoMat.row(j).dot(alphai));
        }
        */
        /*
        VectorXf annoEff = annoMat * alphai;
        #pragma omp parallel for
        for(int j = 0; j < numSnps; j++){
            snpP(j, i) = Normal::cdf_01(annoEff[j]);
        }
        */
    }

    MatrixXf annoEff = annoMat * alpha;
    for (int i = 0; i < numComp; i++) {
        #pragma omp parallel for
        for(int j = 0; j < numSnps; j++){
            snpP(j, i) = Normal::cdf_01(annoEff(j, i));
        }
    }

}

void AnnoProb::compute_AnnoJointProb() {
    computePiFromP(annoCondProb, annoJointProb);
}

void AnnoProb::computePiFromP(const MatrixXf &snpP, MatrixXf &snpPi) {
    int ndist = snpPi.cols();
    int ndist_1 = ndist - 1;

    for (int i = 0; i <= ndist_1; i++) {
        if (i != ndist_1) {
            snpPi.col(i) = 1.0 - snpP.col(i).array();
        } else {
            snpPi.col(i) = VectorXf::Ones(snpP.rows());
        }

        for (int j = 0; j < i; j++) {
            snpPi.col(i) = snpPi.col(i).array() * snpP.col(j).array();
        }
    }
}

VectorXd AnnoProb::getSigmaSq(){
    return sigmaSq.cast<double>();
}

AnnoProb::AnnoProb(string fileAnnot, int numAnno, const VectorXf &Pi, MatrixXf &snpPi, bool bOutDetail, double initSS) {
    this->bOutDetail = bOutDetail;

    numDist = Pi.size();
    numComp = Pi.size() - 1;
    this->numAnno = numAnno;
    numSnps = snpPi.rows();

    // read annot
    annoMat.resize(numSnps, numAnno);
    FILE *hAnnot = fopen(fileAnnot.c_str(), "rb");
    if(!hAnnot){
        Rcout << "Error: can't open the annotation file: " << fileAnnot << std::endl;
        throw("Error");
    }
    int64_t n_ele = (int64_t) numSnps * numAnno;
    if(fread(annoMat.data(), sizeof(float), n_ele, hAnnot) != n_ele){
        Rcout << "Error: read annotation file error, size is not correct: " << fileAnnot << std::endl;
        throw("Error");
    }
    fclose(hAnnot);


    annoSD.resize(numAnno);
    vg_enrich_qt.resize(numAnno); 
    bAnnoBinary.resize(numAnno);

    int numBinary =0;
    int numQT = 0;
    for(int i = 0; i < numAnno; i++){
        bool bBinary = true;
        for(int j = 0; j < numSnps; j++){
            float temp = annoMat(j, i);
            if(abs(temp) > 1e-6 && abs(temp - 1.0) > 1e-6){
                bBinary = false;
                break;
            }
        }
        if(bBinary){
            annoSD[i] = 1.0;
            numBinary++;
        }else{
            annoSD[i] = sqrt((annoMat.col(i).array() - annoMat.col(i).mean()).square().sum() / numSnps);
            numQT++;
        }
        bAnnoBinary[i] = bBinary;
    }
    // get the XPX inverse
    XPXiSD.resize(numAnno);
    XPXqiSD.resize(numAnno);

    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < numAnno; i++){
        float XPX = annoMat.col(i).dot(annoMat.col(i));
        if(XPX < 1e-6){
            Rcout << "Annotation column " << i+1 << " has too small XPX: " << XPX << ", please remove this column from annotation data or perform re-scale." <<std::endl;
            throw("Error");
        }
        if(bAnnoBinary[i]){
            XPXiSD[i] = 1.0 / XPX * annoSD[i];
        }else{
            Matrix2f temp;
            temp(0, 0) = numSnps;
            temp(1, 1) = XPX;
            temp(0, 1) = annoMat.col(i).sum();
            temp(1, 0) = temp(0, 1);
            Eigen::FullPivLU<Matrix2f> lu(temp);
            if(!lu.isInvertible()){
                Rcout << "Annotation column " << i+1 << " can't be inverted, please remove this column from annotation data." <<std::endl;
                throw("Error");
            }
            XPXqiSD[i] = lu.inverse().array() * annoSD[i];
        }
    }
    Rcout << "Total annotation category (including intercept): " << numAnno << "." << std::endl;
    Rcout << "Number of binary annotation: " << numBinary << ", quantitative annotation: " << numQT << "." << std::endl;

    invFrac = numSnps * 1.0 / annoMat.colwise().sum().array();

    if (numSnps != snpPi.rows()) {
        throw("annoMat should have same rows with snpPi");
    }
    if (snpPi.cols() != Pi.size()) {
        throw("snpPi should have same columns with Pi");
    }

    //init the snpP
    snpP.resize(numSnps, numComp);
    // init alpha and snpPi
    alpha = MatrixXf::Zero(numAnno, numComp);
    initP_Pi_annoAlpha(Pi, snpPi);

    // init sigmaSS
    //sigmaSq = VectorXf::Ones(numComp).array() * (float) initSS;
    //sigmaSq = VectorXf::LinSpaced(numComp, 0.1, 1.0) * (float) initSS;
    sigmaSq = VectorXf::Ones(numComp) * (float)initSS;
    /*
    for(int i = 1; i < numComp; i++){
        sigmaSq[i] = sigmaSq[i-1]*2;
    }
    if(sigmaSq[numComp-1] < 1.0) sigmaSq[numComp-1] = 1.0;
    */

    //init ssq
    ssq.resize(numComp);

    annoCondProb.resize(numAnno, numComp);
    annoJointProb.resize(numAnno, numDist);

    annoDists.resize(numAnno, numDist);
    vg_annot_comp.resize(numAnno, numDist);
    vg_tot_annot.resize(numAnno);
    vg_enrich_bin.resize(numAnno);
}

void AnnoProb::samplePi(const MatrixXf &z, MatrixXf &snpPi) {
    //Rcout << "debug: 1" << std::endl;
    annoEffect_sample_Gibbs(z);
    // get snpPi
    //Rcout << "debug: 2" << std::endl;
    computePiFromP(snpP, snpPi);
    snpPi = snpPi.array() + 1e-16;
   // resample sigmaSS
    //Rcout << "debug: 3" << std::endl;
    annoSigmaSS_sample();

}

void AnnoProb::computeProb(){
    // cond prob
    annoCondProb_probit();
    // get joint annot effect
    compute_AnnoJointProb();
}


//compute annote by binary
void AnnoProb::computeEnrichBin(const MatrixXf &vg_snp_comp, float vg){
    //MatrixXf vg_annot_comp(numAnno, numComp);
    //#pragma omp parallel for
    //for(int i = 0; i < numDist; i++){
    //   vg_annot_comp.col(i) = (annoMat.array().colwise() * vg_snp_comp.col(i).array()).colwise().sum();
    //}
    for(int i = 0; i < numDist; i++){
        #pragma omp parallel for
        for(int j = 0; j < numAnno; j++){
            vg_annot_comp(j, i) = annoMat.col(j).dot(vg_snp_comp.col(i));
        }
    }

    vg_tot_annot = vg_annot_comp.rowwise().sum();
    vg_enrich_bin = vg_tot_annot.array() / vg * invFrac.array();
}

// compute enrich by quantitative
void AnnoProb::computeEnrichQt(const VectorXf &snp_vg, float vg){
    if(vg < 1e-5){
        vg_enrich_qt.setZero();
        vg_enrich_qt[0] = 1;
    }else{
        VectorXf XPy = annoMat.transpose() * snp_vg;
        float perSNPh2 = vg / numSnps;
        #pragma omp parallel for schedule(dynamic)
        for(int i = 0; i < numAnno; i++){
            if(bAnnoBinary[i]){
                float coef = XPXiSD[i] * XPy[i];
                vg_enrich_qt[i] = coef / perSNPh2;
            }else{
                VectorXf XPy2(2);
                XPy2(0) = XPy(0);
                XPy2(1) = XPy(i);
                VectorXf coef = XPXqiSD[i] * XPy2;
                vg_enrich_qt[i] = 1.0 + coef(1) / perSNPh2;
            }
        }
    }
}

void AnnoProb::computeDist(const MatrixXf &z, const VectorXf &numSnpMix){
    VectorXi delta = z.rowwise().sum().cast<int>();
    int numDist = numSnpMix.size();
    //MatrixXf annoDists(numAnno, numDist);
    /*
    Rcout << "numAnno: " << annoMat.row(0).size() << ", " << numAnno << std::endl;
    Rcout << "Unique deltas: " << std::endl;
    std::set<int> unique_delta{delta.data(), delta.data() + delta.size()};
    for(auto it = unique_delta.begin(); it != unique_delta.end(); ++it){
        Rcout << ' ' << *it;
    }
    Rcout << std::endl;
    */

    #pragma omp parallel for
    for(int i = 0; i < numDist; i++){
        int nsnpDisti = numSnpMix[i];
        //Rcout << "component " << i << ", nsnpDisti: " << nsnpDisti;
        if(nsnpDisti != 0){
            MatrixXf annotMatCompi(nsnpDisti, numAnno);
            int idx = 0;
            for(int j = 0; j < numSnps; j++){
                if(delta[j] == i){
                    annotMatCompi.row(idx) = annoMat.row(j);
                    idx++;
                }
            }
            //Rcout << ", idx: " << idx << std::endl;

            for(int k = 0; k < numAnno; k++){
                VectorXf annoSrt = annotMatCompi.col(k);
                std::sort(annoSrt.data(), annoSrt.data() + annoSrt.size());
                annoDists(k, i) = annoSrt[annoSrt.size() / 2];
            }

        }else{
            for(int k = 0; k < numAnno; k++){
                annoDists(k, i) = -100000; // very small value; missing
            }
        }
    }
}


