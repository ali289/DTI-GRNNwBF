# DTI-GRNNwBF
Drug-Target Interaction prediction using unifying of graph regularized nuclear norm with bilinear factorization

In this Code, the setting of datasets is based on recent works done on the DTI problem. for example: https://github.com/aanchalMongia/MGRNNMforDTI


Three cross-validation settings (CVS) as named CVS1, CVS2 and CVS3 are implemented. 
In CVS1, which is a common setting for evaluation, the target-drug pairs for the test set were randomly selected for prediction. 
In CVS2 and CVS3, settings are performed to evaluate the ability of methods to predict interactions for novel drugs 
(i.e. drugs for which no interaction information is available) and novel targets, respectively.
It can be pointed out that in CVS2, entire drug profiles and in CVS3, entire target profiles are selected as test set.



Ali Ghanbari sorkhi1, Zahra abbasi2, Majid Iranpour mobarakeh3, Jamshid Pirgazi4

1 Faculty of Electrical and Computer Engineering, University of Science and Technology of Mazandaran, P.O. Box 48518-78195, Behshahr, Iran. ali.ghanbari@mazust.ac.ir
2 Faculty of Medical Biotechnology, School of Medicine, Shahroud University of Medical Sciences, Shahroud, Iran. abbasizahra15@yahoo.com
3 Faculty of Computer Engineering and IT, Payam Noor University. Tehran, Iran, iranpour@pnu.ac.ir
4 Faculty of Electrical and Computer Engineering, University of Science and Technology of Mazandaran, P.O. Box 48518-78195, Behshahr, Iran. j.pirgazi@mazust.ac.ir
