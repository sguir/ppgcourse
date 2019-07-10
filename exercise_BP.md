# Ppgcourse - BayPass


In this practical class we will try to identify human SNPs with evidences of local adaptive diversification and/or significantly associated to geographical variables while taking into account the null correlation of allele frequencies across populations.

* We will use a modified version of the data used in Coop et al. (2010). It consists of genotypes at 2333 + 2 SNPs for 927 individuals from 52 human populations of the Human Genome Diversity Project (HGDP) panel (Conrad et al. 2006). The last two SNPs (rs12913832 and rs1042602) were added manually. 
These two SNPs have been previously reported to be under positive selection in European populations and are located in genes (ERC2 and TYR) involved in light skin and eye color (Wilde et al.2014).
* The covariates (geographical variables) explored are latitude, longitude and a categorical variable with value = 1 if the population is European and -1 if is not.
 

```diff
- Note that even we are performing the analysis with a very reduce dataset, the time it takes to run each model is above what we have in this practice and therefore we will stop BayPass and use some of the results previously obtained to plot them in R and to run the next models.

```

## Get data
The data for this session can be retrieved from google drive [data](https://drive.google.com/file/d/1seD9x50Gf-j7xH7vd7dKMHZOuWfOfMzt/view)   
Download the file "ppg_bp_2019.tar.gz" in the shared folder between the container and the host system (/ppgdata). 
Then, go back to the container terminal and type:

```bash
cd ppgdata
tar -xvzf ppg_bp_2019.tar.gz
```
The folder ppg_bp_2019 has three main subfolders:

* input_data: genotype and covariate input data and the script baypass_utils.R needed for some analysis. 
* results: previously obtained results for a matter of visualization.   
* forR: previously obtained results that are necessary to plot some results obtained during this practical class and to execute some of the BayPass models (since we will not have enough time to run everything during the class).  


The files in each subfolder are classified according to the model/process (e.g., CORE, AUX,...)

Open another container: the one (on your right) will be used to run Baypass ("BayPass container") and the other (on your left) to perform analysis and plots in R (just to avoid to upload the R libraries each time). 

Type in the "BayPass" container terminal:

```bash
cd ppgdata
mkdir baypass

cd ppgdata/ppg_bp_2019/input_data
cp * /ppgdata/baypass/
```

This will copy the input data and the baypass_utils.R script into the folder baypass in which we are going to run baypass

In the "R" container, upload the R libraries:

```R
require(corrplot); require(ape); require(geigen);require(mvtnorm)
source("baypass_utils.R")
```
Now we are ready to run BayPass:

## The CORE Model
The core model allows to perform genome scan for differentiation (covariate free) using the XtX statistics (\~Fst).
The main advantage of this approach is that it explicitly account for the covariance structure in population allele frequencies (via estimating) resulting from the demographic history of the populations.

To run this model with allele data you will need:
* The number of populations in the analysis (```-npop flag```)
* The genotype file (hgdp.geno in PPG_BP_2019/input_data/): the genotypes for each SNP and population. In rows, the SNPs “sorted if possible”. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile flag```). 
* A random number for the seed (in case of needed; ```-seed flag```)
* A prefix to name the output (```-outprefix flag```)

** For running the models with Pool-seq data see the specifications in the manual [BayPass manual](http://www1.montpellier.inra.fr/CBGP/software/baypass/files/BayPass_manual_2.1.pdf)

Run BayPass under the CORE model with three different seeds:

```
g_baypass -npop 52 -gfile hgdp.geno -seed 15263 -outprefix hgdp_s1

g_baypass -npop 52 -gfile hgdp.geno -seed 26847 -outprefix hgdp_s2

g_baypass -npop 52 -gfile hgdp.geno -seed 94875 -outprefix hgdp_s3
```

On the screen, it will apear the specifications of the input file (number of markers, Genotype file name...) and the specifications of the MCMC.

```diff
- Stop BayPass.

```

Copy the previously obtained results to the baypass folder:

```bash
cd ppgdata/ppg_bp_2019/forR/CORE
cp *  /ppgdata/baypass/
```


Upload the estimate of omega (covariance matrix) for each seed:

```R
omega_s1=as.matrix(read.table(file="hgdp_s1_mat_omega.out", header=F))
omega_s2=as.matrix(read.table(file="hgdp_s2_mat_omega.out", header=F))
omega_s3=as.matrix(read.table(file="hgdp_s3_mat_omega.out", header=F))
```

Setup the population names for each omega matrix:

```R
pop.names=c("Papuan","Melanesian","Surui","Pima","Maya","Karitiana","Columbian","Yi","Yaku","Xibo","Uygur","Tujia","Tu","She","Orogen","Naxi","Mongolia","Miao","Lahu","Japanese","Hezhen","Han","Daur","Dal","Cambodian","Sindhi","Pathan","Makrani","Kalash","Hazara","Burusho","Brahui","Balochi","Palestinian","Mozabite","Druze","Bedouin","Tuscan","Sardinian","Russian","Orcadian","French","Italian","Basque","Adygei","Yoruba","San","MbutiPygmy","Mandenka","BiakaPygmy","BantuSouthAfrica","BantuKenya")

dimnames(omega_s1)=list(pop.names,pop.names)
dimnames(omega_s2)=list(pop.names,pop.names)
dimnames(omega_s3)=list(pop.names,pop.names)
```

We can explore the shared history of populations by transforming the covariance matrix into a correlation matrix or a bifurcating phylogenetic tree:

Transform the covariance matrix into a correlation matrix using the R function cov2cor():

```R
cor.mat_s1=cov2cor(omega_s1)
cor.mat_s2=cov2cor(omega_s2)
cor.mat_s3=cov2cor(omega_s3)
```

Plot the correlation matrix (here only for seed1):

```R
pdf(file="correlation_matrix_core_model.pdf")
corrplot(cor.mat_s1, method="color",mar=c(2,1,2,2)+0.1, main=expression("Seed1: Correlation map based on"~hat(Omega)), tl.cex=0.5)
dev.off()
```

Transform the correlation matrix into a hierarchical clustering tree using the R function hclust():

```R
bta14.tree_s1=as.phylo(hclust(as.dist(1-cor.mat_s1**2)))
bta14.tree_s2=as.phylo(hclust(as.dist(1-cor.mat_s2**2)))
bta14.tree_s3=as.phylo(hclust(as.dist(1-cor.mat_s3**2)))
```

Plot the hierarchical clustering tree (here ony for seed1):

```R
pdf(file="correlation_matrix_core_model_tree.pdf")
plot(bta14.tree_s1, type="p",
main=expression("Seed1: Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"), cex=0.5)
dev.off()
```
* We can compare the omega matrices obtained under the CORE model when using different seeds.   
* Remember that we run the CORE model with different seeds to check consistency in the estimation of parameters of the model.

Plot the comparison between omega-seed1 and omega-seed2:
```R
pdf(file="omega_s1_s2_comparison.pdf")
plot(omega_s1, omega_s2) ; abline(a=0,b=1)
dev.off()
```

Compute the distances between pairs of omegas:

```R
dist.12=fmd.dist(omega_s1, omega_s2)
dist.13=fmd.dist(omega_s1, omega_s3)
dist.23=fmd.dist(omega_s2, omega_s3)
```

If there omegas are not significantly different we can assume that there is consistency in the parameters estimation and hence, you should choose one of the omegas to perform the subsequent analyses.


Explore the values of the XtX statistic (~Fst) obtained under the CORE model.

Read the XtX file (here only for seed1):

```R
hgdp_s1.snp.res=read.table("hgdp_s1_summary_pi_xtx.out",h=T)
```

Plot the XtX values:

```R
pdf(file="XtX_core_model_seed1_non-calibrated.pdf")
plot(hgdp_s1.snp.res$M_XtX, xlab="SNPs", ylab="XtX", main="Seed 1")
dev.off()
```
```QUESTION: which XtX values are significantly different?```

### Pseudo Observed Data (PODs) 
Here, we are going to simulate data (PODs) using the R function simulate.baypass() in the baypass_utils.R script (provided in the BayPass package).
PODs are simulated under the inference model (e.g., using posterior estimates of the covariance matrix and the a and b parameters of the beta prior distribution for the
overall (across population) SNP allele frequencies).
Once these PODS are simulated, we need to run again the CORE model to built the \"expected\" distribution of the XtX values under the inference model in order to find which of the observed XtX values are significantly different from the expected (calibration process)


Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution:

```R
pi.beta.coef=read.table("hgdp_s1_summary_beta_params.out",h=T)$Mean
```

Upload the original data to obtain total allele count (sample size for each population). 
Do this by using the geno2YN() function in baypass_utils.R script:

```R
hgdp.data<-geno2YN("hgdp.geno")
```

Read the omega matrix from seed1:

```R
omega_s1=as.matrix(read.table(file="hgdp_s1_mat_omega.out", header=F))
```


Run **only the first** (simu.hgdp_1000) of the two simulation experiments with different number of replicates (experiment1=1000 and experiment2=100000 replicates)  
* We are not going to run the second one for a matter of time
```R
simu.hgdp_1000 <- simulate.baypass(omega.mat=omega_s1, nsnp=1000, sample.size=hgdp.data$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="hgdp_pods_1000")

simu.hgdp_100000 <- simulate.baypass(omega.mat=omega_s1, nsnp=100000, sample.size= hgdp.data$NN, beta.pi=pi.beta.coef, pi.maf=0, suffix="hgdp_pods_100000")
```

* Note that the G.hgdp_pods_1000 and G.hgdp_pods_100000 files are now the new genotype input files resulting from the simulation process.  


Run again the Core model **only with the first simulation experiment** (G.hgdp_pods_1000) (remember, in the BayPass container!)

```
g_baypass -npop 52 -gfile G.hgdp_pods_1000 -outprefix hgdp_pod_1000 

g_baypass -npop 52 -gfile G.hgdp_pods_100000 -outprefix hgdp_pod_100000

```

```diff
- Stop BayPass.

```

Copy the previously obtained results to the baypass folder:

```bash
cd ppgdata/ppg_bp_2019/forR/PODS/CORE/1000
cp *  /ppgdata/baypass/

cd ppgdata/ppg_bp_2019/forR/PODS/CORE/100000
cp *  /ppgdata/baypass/

```

### Sanity Check 
Here, we are comparing the simulated data (PODS) under the inference model to the observed data to assess if the inference model (posteriors for the covariance matrix and the other hyperparameters) are giving us \"valid\" predictions about the \"reality\".
In other words, if the model we have inferred is able to generate data similar to the observed ones.


#### 1000 PODS
Get the omega estimated from the POD analysis:
```R
pod.omega_1000=as.matrix(read.table("hgdp_pod_1000_mat_omega.out"))
```

Plot the observed versus the simulated omegas:
```R
pdf(file="comparison_omega_obs_sim_1000pods.pdf")
plot(pod.omega_1000, omega_s1) ; abline(a=0,b=1)
dev.off()
```

Compute the distance between the observed and simulated omegas
```R
fmd.dist(pod.omega_1000, omega_s1)
```

Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the PODs:
```R
pod.pi.beta.coef_1000=read.table("hgdp_pod_1000_summary_beta_params.out",h=T)$Mean
```

Plot the observed versus the simulated pi.beta:
```R
pdf(file="comparison_pi.beta_obs_sim_1000pods.pdf")
plot(pod.pi.beta.coef_1000, pi.beta.coef) ; abline(a=0,b=1)
dev.off()
```

#### 100000 PODS
Get the omega estimated from the POD analysis:
```R
pod.omega_100000=as.matrix(read.table("hgdp_pod_100000_mat_omega.out"))
```

Plot the observed versus the simulated omegas:
```R
pdf(file="comparison_omega_obs_sim_100000pods.pdf")
plot(pod.omega_100000, omega_s1) ; abline(a=0,b=1)
dev.off()
```

Distance between the observed and the simulated omegas:
```R
fmd.dist(pod.omega_100000, omega_s1)
```

Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the PODs:
```R
pod.pi.beta.coef_100000=read.table("hgdp_pod_100000_summary_beta_params.out",h=T)$Mean
```

Plotting the observed versus the simulated pi.beta:
```R
pdf(file="comparison_pi.beta_obs_sim_100000pods.pdf")
plot(pod.pi.beta.coef_100000, pi.beta.coef) ; abline(a=0,b=1)
dev.off()
```
```QUESTION: What is the main difference when comparing the two simulation experiments (1000 and 1000000 PODs) to the observed data?```

### XtX calibration 
Remember that this is for building the \"expected\" distribution of the XtX values under the inference model in order to find which of the observed XtX values are significantly different from the expected 

Get the POD XtX's (1000 pods):
```R
pod.xtx=read.table("hgdp_pod_1000_summary_pi_xtx.out", h=T)$M_XtX
```
Compute the 1% threshold using the R function quantile(): 
```R
pod.thresh=quantile(pod.xtx, probs=0.99)
```

Add the threshold to the XtX observed data plot: 
```R
pdf(file="XtX_calibration_core_model_1000.pdf")
plot(hgdp_s1.snp.res$M_XtX, xlab="SNPs", ylab=" XtX ")
abline(h=pod.thresh, lty=2, col="red")
dev.off()
```

Retrive the SNPs that are significant (XtX): 
```R
snps.sign.xtx <- hgdp_s1.snp.res[hgdp_s1.snp.res$M_XtX >= pod.thresh, ]
write.table(snps.sign.xtx, file="significant_SNPs_XtX.txt", col.names=T, row.names=F, quote=F, sep="\t")
```
```QUESTION: Which are the SNPs that are significant?```  
Here [map file]() you can find a "map" file associated to you original geno file with the name of each SNP (chromosome and position can be found in NCBI database)

## The AUX Model
This model allows to evaluate to which extent the population covariables are (linearly) associated to each marker/SNP.
The introduction of an auxiliary variable (delta) in the AUX model indicates whether a specific SNP i can
be regarded as associated to a given covariable k (delta_ik = 1) or not (delta_ik = 0).
We can compare both models (association and no-association) through Bayes Factor (BFmc) while dealing with multiple testing issues.

It is (strongly) recommended to scale each covariate (so that the mean = 0 and the std = 1 for each covariable). 
The scalecov option allows to perform this step automatically prior to analysis, if needed.

To run this model with allele data you will need:
* The number of populations in the analysis (```-npop flag```)
* The genotype file (hgdp.geno in PPG_BP_2019/input_data/): the genotypes for each SNP and population. In rows, the SNPs “sorted if possible”. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile flag```).
* The covariates file (covariates in PPG_BP_2019/input_data/): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile flag```).
* To specify if you want to scale covariables (```-scalecov flag```)
* To specify the model you want to run (```-auxmodel flag```)
* To provide the omega matrix obtained under the CORE model (```-omegafile flag```)
* A prefix to name the output (```-outprefix flag```)

Run the AUX model (remember, in the BayPass container!):

```
g_baypass -npop 52 -gfile hgdp.geno -efile covariates -scalecov -auxmodel -omegafile hgdp_s1_mat_omega.out -outprefix hgdpaux
```

On the screen, it will apear the specifications of the input file (number of markers, Genotype file name, Covariables file...) and the specifications of the MCMC.


```diff
- Stop BayPass.

```

Copy the previously obtained results to the baypass folder:

```bash
cd ppgdata/ppg_bp_2019/forR/AUX/
cp *  /ppgdata/baypass/

```

Get the regression coefficients (posterior mean) for each SNP:
```R
covaux.snp.res=read.table("hgdpaux_summary_betai.out",h=T)
```

Plot the resulting estimates of the Bayes Factor, the underlying regression coefficients (posterior mean) and the calibrated XtX:

```R
pdf(file="plot_Aux_model_results.pdf")
layout(matrix(1:3,3,1))
plot(covaux.snp.res$BF.dB., xlab="SNP",ylab="BFmc (in dB)")
abline(h=10, lty=2, col="green")
abline(h=15, lty=2, col="blue")
abline(h=20, lty=2, col="red")
plot(covaux.snp.res$M_Beta, xlab="SNP",ylab=expression(beta~"coefficient"))
plot(hgdp_s1.snp.res$M_XtX, xlab="SNP",ylab="XtX corrected for SMS")
abline(h=pod.thresh, lty=2, col="red")
dev.off()
```
* Blue, green and red lines in the BF plot are depicting the thresholds for the SNPs with strong, very strong and decisive evidence of association, respectively.
* Note that the calibrated XtX values were obtained from the CORE model and the PODs.
* Note also that the X axis in the two first plots ranges from 0 to more than 7000 SNPs. This is because you are plotting the 2335 SNPs for each of the three covariables.

Retrieve the SNPs with strong (se), very strong (vse) and decisive evidence (de) of significant association:
```R
se.bf <- covaux.snp.res[covaux.snp.res$BF.dB. > 10 & covaux.snp.res$BF.dB. < 15, ]
vse.bf <- covaux.snp.res[covaux.snp.res$BF.dB. > 15 & covaux.snp.res$BF.dB. < 20, ]
de.bf <- covaux.snp.res[covaux.snp.res$BF.dB. >= 20, ]

write.table(se.bf, file="aux_model_strong_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(vse.bf, file="aux_model_very_strong_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(de.bf, file="aux_model_decisive_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
```
```QUESTION: Which are the significant SNPs and to which covariable they are associated? Are some of the significant SNPs associated to more than one covariable?```  

## The AUX Model with Linkage Disequilibrium (LD)
This model is similar to the AUX model but in case you know the physical order of the SNPs in chromosomes it allows to account for
spatial dependency among markers (\~linkage disequilibrium) by means of the parameter b_is.
b_is = 0 implies no spatial dependency and b_is > 0 assumes that the delta_ik with similar values tend to cluster according to the underlying SNP positions (the higher the b_is, the higher the level of
spatial homogeneity). In practice, b_is = 1 is commonly used and a value of b_is =< 1 is recommended.

To run this model with allele data you will need:
* The number of populations in the analysis (```-npop flag```)
* The genotype file (hgdp.geno in PPG_BP_2019/input_data/): the genotypes for each SNP and population. In rows, the SNPs “sorted if possible”. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile flag```).
* The covariates file (covariates in PPG_BP_2019/input_data/): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile flag```).
* To set the value of the  b_is parameter (```-isingbeta flag```)
* To specify if you want to scale covariables (```-scalecov flag```)
* To specify the model you want to run (```-auxmodel flag```)
* To provide the omega matrix obtained under the CORE model (```-omegafile flag```)
* A prefix to name the output (```-outprefix flag```)

Run Baypass with the Aux model (in the BayPass container!):
```
g_baypass -npop 52 -gfile hgdp.geno -efile covariates -auxmodel -isingbeta 1.0 -omegafile hgdp_s1_mat_omega.out -outprefix hgdpauxld
```

```diff
- Stop BayPass.

```

Copy the previously obtained results to the baypass folder:

```bash
cd ppgdata/ppg_bp_2019/forR/AUX_LD/
cp *  /ppgdata/baypass/

```

Get the file with the regression coefficients (posterior mean) for each SNP: 
```R
covauxisb1.snp.res=read.table("hgdpauxld_summary_betai.out",h=T)
```

Plot the estimates of the posterior mean of the each auxiliary variable delta_ik under both models (Aux and Aux with LD):

```R
pdf(file="plot_Aux_ld_model_results_delta.pdf")
layout(matrix(1:2,2,1))
plot(covaux.snp.res$M_Delta,xlab="SNP",ylab=expression(delta[i]),main="AUX model")
plot(covauxisb1.snp.res$M_Delta,xlab="SNP",ylab=expression(delta[i]), main="AUX model with isb=1")
dev.off()
```
* Note also that the X axis in the two first plots ranges from 0 to more than 7000 SNPs. This is because you are plotting the 2335 SNPs for each of the three covariables.

Plot the resulting estimates of the Bayes Factor, the underlying regression coefficients (posterior mean) and the calibrated XtX:

```R
pdf(file="plot_Aux_ld_model_results.pdf")
layout(matrix(1:3,3,1))
plot(covauxisb1.snp.res$BF.dB., xlab="SNP",ylab="BFmc (in dB)")
abline(h=10, lty=2, col="green")
abline(h=15, lty=2, col="blue")
abline(h=20, lty=2, col="red")
plot(covauxisb1.snp.res$M_Beta, xlab="SNP",ylab=expression(beta~"coefficient"))
plot(hgdp_s1.snp.res$M_XtX, xlab="SNP",ylab="XtX corrected for SMS")
abline(h=pod.thresh, lty=2, col="red")
dev.off()
```

* Again, blue, green and red lines in the BF plot are depicting the thresholds for the SNPs with strong, very strong and decisive evidence of association, respectively.
* Note that the calibrated XtX values were obtained from the CORE model and the PODs.
* Note also that the X axis in the two first plots ranges from 0 to more than 7000 SNPs. This is because you are plotting the 2335 SNPs for each of the three covariables.

Retrieve the SNPs with strong (se), very strong (vse) and decisive evidence (de) of significant association:
```R
se.bf <- covauxisb1.snp.res[covauxisb1.snp.res$BF.dB. > 10 & covauxisb1.snp.res$BF.dB. < 15, ]
vse.bf <- covauxisb1.snp.res[covauxisb1.snp.res$BF.dB. > 15 & covauxisb1.snp.res$BF.dB. < 20, ]
de.bf <- covauxisb1.snp.res[covauxisb1.snp.res$BF.dB. >= 20, ]

write.table(se.bf, file="aux_ld_model_strong_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(vse.bf, file="aux_ld_model_very_strong_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(de.bf, file="aux_ld_model_decisive_evidence_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
```

```QUESTION: Which are the SNPs that are significantly associated? is the number of significant SNPs more or less the same as that obtained under the AUX model?```

# EXTRA EXERCISES

## The STANDARD Model: importance sampling 
This model allows to evaluate to which extent the population covariables are (linearly) associated to each marker/SNP.
The estimation of the beta regression coefficients for each SNP and covariable is performed using the importance sampling approach.

* Remember that this model is recommended when the number of populations is small (e.g.,  8) and/or when populations are highly differentiated.
* The importance sampling algorithm relies on approximations to estimate the regression coefficients (do not sample them from the posterior distribution and hence, the model should be run 3-5 times with different seeds to check consistency across runs). 
* The median of beta among the runs can be taken as the estimate beta parameter.
* Bayes factor importance sampling (BFis) and the approximated Bayesian P-value (eBPis) are used to evaluate if a particular SNP is associated with a particular covariable and as the XtX statistic in the CORE model, they should be calibrated. 
* Note that unlike the other models (AUX, AUX-LD, STDmcmc), this model calculates the covariance matrix (omega) and the correlation parameter (beta) at the same time.

To run this model with allele data you will need:
* The number of populations in the analysis (```-npop flag```)
* The genotype file (hgdp.geno in PPG_BP_2019/input_data/): the genotypes for each SNP and population. In rows, the SNPs “sorted if possible”. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile flag```).
* The covariates file (covariates in PPG_BP_2019/input_data/): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile flag```).
* To specify if you want to scale covariables (```-scalecov flag```)
* A prefix to name the output (```-outprefix flag```)

Run BayPass with the STANDARD model importance sampling: 

```
g_baypass -npop 52 -gfile hgdp.geno -efile covariates -scalecov -outprefix hgdpis
```

```diff
- Stop BayPass.

```

Copy the previously obtained results to the baypass folder:

```bash
cd ppgdata/ppg_bp_2019/forR/STDis/
cp *  /ppgdata/baypass/

```

### Pseudo Observed Data (PODs) 
Here and for a matter of time, we are going to simulate 1,000 PODs (but usually 100000 is recommended). 
To simulate PODS under the STDis we need to follow the same procedure as with XtX in the CORE model but adding to extra parameters: a vector with the values of EACH covariable and a vector with the estimated mean beta parameter for EACH covariable (output resulting from running the STDis model). 

Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution for each covariable separately: 

```R
pi.beta.coefis=read.table("hgdpis_summary_beta_params.out",h=T)$Mean
```

Get estimates (posterior mean) of the correlation parameter (beta) for each covariable separately:

```R
beta.coefis=read.table("hgdpis_summary_betai_reg.out",h=T)
beta.coefislat=beta.coefis[beta.coefis$COVARIABLE == 1, ]$Beta_is
beta.coefislon=beta.coefis[beta.coefis$COVARIABLE == 2, ]$Beta_is 
beta.coefiseu=beta.coefis[beta.coefis$COVARIABLE == 3, ]$Beta_is
```

Upload the original data to obtain total allele count (sample size for each population):

```R
hgdp.data<-geno2YN("hgdp.geno")
```
Get the omega matrix estimated under the STDis model:

```R
omega.is=as.matrix(read.table(file="hgdpis_mat_omega.out", header=F))
```

Get the values for each covariable separately:

```R
cov=read.table("covariates",h=F)
cov1=as.numeric(cov[1,])
cov2=as.numeric(cov[2,])
cov3=as.numeric(cov[3,])
```

Create the PODs for each covariable separately:
```R
simu.hgdpislat_1000 <- simulate.baypass(omega.mat=omega.is, nsnp=1000, beta.coef=beta.coefislat, beta.pi=pi.beta.coefis, pop.trait=cov1, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsislat_1000")

simu.hgdpislon_1000 <- simulate.baypass(omega.mat=omega.is, nsnp=1000, beta.coef=beta.coefislon, beta.pi=pi.beta.coefis, pop.trait=cov2, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsislon_1000")

simu.hgdpiseu_1000 <- simulate.baypass(omega.mat=omega.is, nsnp=1000, beta.coef=beta.coefiseu, beta.pi=pi.beta.coefis, pop.trait=cov3, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsiseu_1000")
```

Run again the STDis model with the simulated PODs and for each covariable independently:

```
g_baypass -npop 52 -gfile G.hgdp_podsislat_1000 -efile covariates_lat -scalecov -outprefix pod_hgdpvislat_1000

g_baypass -npop 52 -gfile G.hgdp_podsislon_1000 -efile covariates_lon -scalecov -outprefix pod_hgdpvislon_1000

g_baypass -npop 52 -gfile G.hgdp_podsiseu_1000 -efile covariates_eu -scalecov -outprefix pod_hgdpviseu_1000
```

```diff
- Stop BayPass.

```

```bash
cd ppgdata/ppg_bp_2019/forR/PODS/STDis/
cp *  /ppgdata/baypass/

```

* Note that for running the STDis a file with all the covariables is normally used whereas to run the simulations under the STDis each variable should be included at time (Number of runs of simulations = Number of covariables). 

### XtX, BFis and eBPis calibration  

Get the pod XtX (1000 pods):

```R
podis.xtx.lat=read.table("pod_hgdpvislat_1000_summary_pi_xtx.out", h=T)$M_XtX
hgdp.is.snp.res.lat <- read.table("hgdpis_summary_pi_xtx.out", h=T)

podis.xtx.lon=read.table("pod_hgdpvislon_1000_summary_pi_xtx.out", h=T)$M_XtX
hgdp.is.snp.res.lon <- read.table("hgdpis_summary_pi_xtx.out", h=T)

podis.xtx.eu=read.table("pod_hgdpviseu_1000_summary_pi_xtx.out", h=T)$M_XtX
hgdp.is.snp.res.eu <- read.table("hgdpis_summary_pi_xtx.out", h=T)
```

Compute the 1% threshold:

```R
podis.thresh.lat=quantile(podis.xtx.lat, probs=0.99)
podis.thresh.lon=quantile(podis.xtx.lon, probs=0.99)
podis.thresh.eu=quantile(podis.xtx.eu, probs=0.99)
```

Add the threshold to the actual XtX plot:

```R
pdf(file="XtX_calibration_stdis.lat.pdf")
plot(hgdp.is.snp.res.lat$M_XtX, xlab="SNPs", ylab=" XtX ")
abline(h=podis.thresh.lat, lty=2, col="red")
dev.off()

pdf(file="XtX_calibration_stdis.lon.pdf")
plot(hgdp.is.snp.res.lon$M_XtX, xlab="SNPs", ylab=" XtX ")
abline(h=podis.thresh.lon, lty=2, col="red")
dev.off()

pdf(file="XtX_calibration_stdis.eu.pdf")
plot(hgdp.is.snp.res.eu$M_XtX, xlab="SNPs", ylab=" XtX ")
abline(h=podis.thresh.eu, lty=2, col="red")
dev.off()
```

Retrive the SNPs that are significant (XtX statsitic):

```R
snps.sign.xtx.is.lat <- hgdp.is.snp.res.lat[hgdp.is.snp.res.lat$M_XtX >= podis.thresh.lat, ]
write.table(snps.sign.xtx.is.lat, file="significant_SNPs_XtX.std.is.lat.txt", col.names=T, row.names=F, quote=F, sep="\t")

snps.sign.xtx.is.lon <- hgdp.is.snp.res.lon[hgdp.is.snp.res.lon$M_XtX >= podis.thresh.lon, ]
write.table(snps.sign.xtx.is.lon, file="significant_SNPs_XtX.std.is.lon.txt", col.names=T, row.names=F, quote=F, sep="\t")

snps.sign.xtx.is.eu <- hgdp.is.snp.res.eu[hgdp.is.snp.res.eu$M_XtX >= podis.thresh.eu, ]
write.table(snps.sign.xtx.is.eu, file="significant_SNPs_XtX.std.is.eu.txt", col.names=T, row.names=F, quote=F, sep="\t")
```

Get the resulting regression coefficient beta, the Importance Sampling estimates of the Bayes Factor, the empirical Bayesian P-value and the underlying regression coefficient for each covariable:

```R
covis.snp.res=read.table("hgdpis_summary_betai_reg.out",h=T)
```

Compute the 1% threshold for BFis and eBPis and for each covariable:

```R

podis.ebp.lat=read.table("pod_hgdpvislat_1000_summary_betai_reg.out",h=T)
podis.thresh.ebp.lat=quantile(podis.ebp.lat$eBPis, probs=0.99)
bfis.lat=podis.ebp.lat$BF.dB
bfis.thresh.lat=quantile(bfis.lat, probs=0.99)

podis.ebp.lon=read.table("pod_hgdpvislon_1000_summary_betai_reg.out",h=T)
podis.thresh.ebp.lon=quantile(podis.ebp.lon$eBPis, probs=0.99)
bfis.lon=podis.ebp.lon$BF.dB
bfis.thresh.lon=quantile(bfis.lon, probs=0.99)

podis.ebp.eu=read.table("pod_hgdpviseu_1000_summary_betai_reg.out",h=T)
podis.thresh.ebp.eu=quantile(podis.ebp.eu$eBPis, probs=0.99)
bfis.eu=podis.ebp.eu$BF.dB.
bfis.thresh.eu=quantile(bfis.eu, probs=0.99)
```

Plot the resulting Importance Sampling estimates of the Bayes Factor, the empirical Bayesian P-value and the underlying regression coefficient for each covariable:

```R
pdf(file="plot_stdis_model_results_lat.pdf")
layout(matrix(1:3,3,1))
plot(covis.snp.res[covis.snp.res$COVARIABLE == 1, ]$BF.dB., xlab="SNP",ylab="BFis (in dB)")
abline(h=bfis.thresh.lat, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 1, ]$eBPis, xlab="SNP",ylab="eBPis")
abline(h=podis.thresh.ebp.lat, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 1, ]$Beta_is, xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

pdf(file="plot_stdis_model_results_lon.pdf")
layout(matrix(1:3,3,1))
plot(covis.snp.res[covis.snp.res$COVARIABLE == 2, ]$BF.dB., xlab="SNP",ylab="BFis (in dB)")
abline(h=bfis.thresh.lon, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 2, ]$eBPis, xlab="SNP",ylab="eBPis")
abline(h= podis.thresh.ebp.lon, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 2 ,]$Beta_is, xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()

pdf(file="plot_stdis_model_results_eu.pdf")
layout(matrix(1:3,3,1))
plot(covis.snp.res[covis.snp.res$COVARIABLE == 3, ]$BF.dB., xlab="SNP",ylab="BFis (in dB)")
abline(h=bfis.thresh.eu, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 3, ]$eBPis, xlab="SNP",ylab="eBPis")
abline(h= podis.thresh.ebp.eu, lty=2, col="red")
plot(covis.snp.res[covis.snp.res$COVARIABLE == 3 ,]$Beta_is, xlab="SNP",ylab=expression(beta~"coefficient"))
dev.off()
```

Retrieve the SNPs with significant association (BFis, eBP) for each covariable:

```R
bf.is.lat <- covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.lat, ][covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.lat, ]$COVARIABLE == 1, ]
ebp.is.lat <- covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.lat, ][covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.lat, ]$COVARIABLE == 1, ]

bf.is.lon <- covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.lon, ][covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.lon, ]$COVARIABLE == 2, ]
ebp.is.lon <- covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.lon, ][covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.lon, ]$COVARIABLE == 2, ]

bf.is.eu <- covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.eu, ][covis.snp.res[covis.snp.res$BF.dB. > bfis.thresh.eu, ]$COVARIABLE == 3, ]
ebp.is.eu <- covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.eu, ][covis.snp.res[covis.snp.res$eBPis > podis.thresh.ebp.eu, ]$COVARIABLE == 3, ]

write.table(bf.is.lat, file="std.is_signif_BFis_lat_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(ebp.is.lat, file="std.is_signif_eBP_lat_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")

write.table(bf.is.lon, file="std.is_signif_BFis_lon_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(ebp.is.lon, file="std.is_signif_eBP_lon_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")

write.table(bf.is.eu, file="std.is_signif_BFis_eu_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(ebp.is.eu, file="std.is_signif_eBP_eu_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
```
## The STANDARD Model: mcmc
This model allows to evaluate to which extent the population covariables are (linearly) associated to each marker/SNP.
The estimation of the beta regression coefficients for each SNP and covariable is performed using an mcmc approach.

* In this case, the user should provide the omega matrix (e.g., using posterior estimates obtained from the CORE model) and it is recommended to
consider only one covariable at a time (particularly if some covariables are correlated).
* Bayes factor mcmc (eBPmcmc) is used to evaluate if a particular SNP is associated with a particular covariable. And as the XtX statistic in the CORE model, they should be calibrated. 

To run this model with allele data you will need:
* The number of populations in the analysis (```-npop flag```)
* The genotype file (hgdp.geno in PPG_BP_2019/input_data/): the genotypes for each SNP and population. In rows, the SNPs “sorted if possible”. In columns: populations. Each population has two columns: one for the reference and the other for the alternative allele counts (```-gfile flag```).
* The covariates file (covariates in PPG_BP_2019/input_data/): In rows, the covariates. In columns, populations (one column per population). The order of the populations should be the same as in the genotype file (```-efile flag```).
* To specify the model you want to run (```-covmcmc flag```)
* To specify if you want to scale covariables (```-scalecov flag```)
* To provide the omega matrix obtained under the CORE model (```-omegafile flag```)
* A prefix to name the output (```-outprefix flag```)

Run STANDARD model independently with each covariable:

```
g_baypass -npop 52 -gfile hgdp.geno -efile covariates_lat -covmcmc -scalecov -omegafile hgdp_s1_mat_omega.out -outprefix hgdp_mcmc_lat

g_baypass -npop 52 -gfile hgdp.geno -efile covariates_lon -covmcmc -scalecov -omegafile hgdp_s1_mat_omega.out -outprefix hgdp_mcmc_lon

g_baypass -npop 52 -gfile hgdp.geno -efile covariates_eu -covmcmc -scalecov -omegafile hgdp_s1_mat_omega.out -outprefix hgdp_mcmc_eu
```

```diff
- Stop BayPass. 
```


```bash
cd ppgdata/ppg_bp_2019/forR/STDmcmc/
cp *  /ppgdata/baypass/

```

### Pseudo Observed Data (PODs) 
Again, for a matter of time, we are going to simulate 1,000 PODs (but usually 100000 is recommended). 
To simulate PODS under the STDmcmc we need to follow the same procedure as with XtX in the CORE model but adding to extra parameters: a vector with the values of EACH covariable and a vector with the estimated mean beta parameter for EACH covariable (output resulting from running the STDmcmc model). 


Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution:

```R
pi.beta.coefmc.lat=read.table("hgdp_mcmc_lat_summary_beta_params.out", h=T)$Mean
pi.beta.coefmc.lon=read.table("hgdp_mcmc_lon_summary_beta_params.out", h=T)$Mean
pi.beta.coefmc.eu=read.table("hgdp_mcmc_eu_summary_beta_params.out", h=T)$Mean
```

Upload the original data to obtain total allele count (sample size for each population)

```
hgdp.data<-geno2YN("hgdp.geno")
```

Get the omega matrix estimated under the CORE model:

```R
omega.mc.lat=as.matrix(read.table(file="hgdp_s1_mat_omega.out", header=F))
```

Get estimates (posterior mean) of the correlation parameter (beta) for each covariable separately: 

```R
beta.coefmc.lat=read.table("hgdp_mcmc_lat_summary_betai.out",h=T)$M_Beta
beta.coefmc.lon=read.table("hgdp_mcmc_lon_summary_betai.out",h=T)$M_Beta
beta.coefmc.eu=read.table("hgdp_mcmc_eu_summary_betai.out",h=T)$M_Beta
```

Get the values for each covariable separately:
cov=read.table("covariates",h=F)
cov1=as.numeric(cov[1,])
cov2=as.numeric(cov[2,])
cov3=as.numeric(cov[3,])

Create the PODs for each covariable separately:

```R
simu.hgdpmclat_1000 <- simulate.baypass(omega.mat=omega.mc.lat, nsnp=1000, beta.coef=beta.coefmc.lat, beta.pi=pi.beta.coefmc.lat, pop.trait=cov1, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsmclat_1000")

simu.hgdpmclon_1000 <- simulate.baypass(omega.mat=omega.mc.lat, nsnp=1000, beta.coef=beta.coefmc.lon, beta.pi= pi.beta.coefmc.lon, pop.trait=cov2, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsmclon_1000")

simu.hgdpmceu_1000 <- simulate.baypass(omega.mat=omega.mc.lat, nsnp=1000, beta.coef=beta.coefmc.eu, beta.pi= pi.beta.coefmc.eu, pop.trait=cov3, sample.size=hgdp.data$NN, pi.maf=0, suffix="hgdp_podsmceu_1000")
```

Run again the STDmcmc model with the simulated PODs and for each covariable independently:

```
g_baypass -npop 52 -gfile G.hgdp_podsmclat_1000 -efile covariates_lat -scalecov -covmcmc -omegafile hgdp_s1_mat_omega.out -outprefix pod_hgdp_mcmc_lat

g_baypass -npop 52 -gfile G.hgdp_podsmclon_1000 -efile covariates_lon -scalecov -covmcmc -omegafile hgdp_s1_mat_omega.out -outprefix pod_hgdp_mcmc_lon

g_baypass -npop 52 -gfile G.hgdp_podsmceu_1000 -efile covariates_eu -scalecov -covmcmc -omegafile hgdp_s1_mat_omega.out -outprefix pod_hgdp_mcmc_eu
```

```diff
- Stop BayPass. 

```

```bash
cd ppgdata/ppg_bp_2019/forR/PODS/STDmcmc/
cp *  /ppgdata/baypass/

```

Get the resulting regression coefficient beta and the Bayes Factor mcmcm (BFmcmc):

```R
covmcmc.snp.res.lat=read.table("hgdp_mcmc_lat_summary_betai.out",h=T)
covmcmc.snp.res.lon=read.table("hgdp_mcmc_lon_summary_betai.out",h=T)
covmcmc.snp.res.eu=read.table("hgdp_mcmc_eu_summary_betai.out",h=T)
```

Compute the 1% threshold for Bayes Factor mcmcm (BFmcmc) for each covariable:

```R
podmc.ebp.lat=read.table("pod_hgdp_mcmc_lat_summary_betai.out",h=T)
podmc.thresh.ebp.lat=quantile(podmc.ebp.lat$eBPmc, probs=0.99)

podmc.ebp.lon=read.table("pod_hgdp_mcmc_lon_summary_betai.out",h=T)
podmc.thresh.ebp.lon=quantile(podmc.ebp.lon$eBPmc, probs=0.99)

podmc.ebp.eu=read.table("pod_hgdp_mcmc_eu_summary_betai.out",h=T)
podmc.thresh.ebp.eu=quantile(podmc.ebp.eu$eBPmc, probs=0.99)
```

Plot the resulting estimates of the empirical Bayesian P-values, the underlying regression coefficients (posterior mean) and the calibrated XtX obtained under the CORE model for each covariable:

```R
pdf(file="plot_std_mcmc_model_results_lat.pdf")
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res.lat$eBPmc, xlab="SNP",ylab="eBPmc")
abline(h=podmc.thresh.ebp.lat, lty=2, col="red")
plot(covmcmc.snp.res.lat$M_Beta, xlab="SNP",ylab=expression(beta~"coefficient"))
plot(hgdp_s1.snp.res$M_XtX, xlab="SNP",ylab="XtX corrected for SMS") 
abline(h=pod.thresh, lty=2, col="red")
dev.off()

pdf(file="plot_std_mcmc_model_results_lon.pdf")
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res.lon$eBPmc, xlab="SNP",ylab="eBPmc")
abline(h=podmc.thresh.ebp.lon, lty=2, col="red")
plot(covmcmc.snp.res.lon$M_Beta, xlab="SNP",ylab=expression(beta~"coefficient"))
plot(hgdp_s1.snp.res$M_XtX, xlab="SNP",ylab="XtX corrected for SMS") 
abline(h=pod.thresh, lty=2, col="red")
dev.off()

pdf(file="plot_std_mcmc_model_results_eu.pdf")
layout(matrix(1:3,3,1))
plot(covmcmc.snp.res.eu$eBPmc, xlab="SNP",ylab="eBPmc")
abline(h=podmc.thresh.ebp.eu, lty=2, col="red")
plot(covmcmc.snp.res.eu$M_Beta, xlab="SNP",ylab=expression(beta~"coefficient"))
plot(hgdp_s1.snp.res$M_XtX, xlab="SNP",ylab="XtX corrected for SMS") 
abline(h=pod.thresh, lty=2, col="red")
dev.off()
```

Retrieve the SNPs with significant association (eBPmc)

```R
bf.mc.lat <- covmcmc.snp.res.lat[covmcmc.snp.res.lat$eBPmc > podmc.thresh.ebp.lat, ]
bf.mc.lon <- covmcmc.snp.res.lon[covmcmc.snp.res.lon$eBPmc > podmc.thresh.ebp.lon, ]
bf.mc.eu <- covmcmc.snp.res.eu[covmcmc.snp.res.eu$eBPmc > podmc.thresh.ebp.eu, ]

write.table(bf.mc.lat, file="std.mcmc_signif_eBPmc_lat_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(bf.mc.lon, file="std.mcmc_signif_eBPmc_lon_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
write.table(bf.mc.eu, file="std.mcmc_signif_eBPmc_eu_SNPs.txt", col.names=T, row.names=F, quote=F, sep= "\t")
```
### BIBLIOGRAPHY

* Gautier M. 2015. Genome-Wide Scan for Adaptive Divergence and Association with Population-Specific Covariates. Genetics 201(4): 1555–1579. 
* Günther T and Coop G. 2013. Robust Identification of Local Adaptation from Allele Frequencies. Genetics 195(1): 205–220.
* Coop G. Witonsky D, Di Rienzo A, Pritchard JK. 2010. Using environmental correlations to identify loci underlying local adaptation. Genetics 185: 1411–1423.
* Hoban S, Kelley JL, Lotterhos KE, Antolin MF, Bradburd G, Lowry DB, Poss ML, Reed LK, Storfer A, Whitlock MC. 2016. Finding the Genomic Basis of Local Adaptation: Pitfalls, Practical Solutions, and Future Directions. Am Nat. 188(4): 379–397



