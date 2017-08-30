# GaussianSV

INTRODUCTION

This is an implementation of GaussianSV - a semantic visualization method for short texts from Le & Lauw (IJCAI 2017).

Usage:

	perl gaussiansv.pl	--num_topics $num_topics
				--dim $dim
				--beta $beta
				--gamma $gamma
				--EM_iter $EM_iter
				--Quasi_iter $Quasi_iter
				--data $data
				--word_vectors $word_vectors_file
				--word_vectors_dim $D
				--rho $rho
				--rho_0 $rho_0
				--output_file $output_file

Arguments:
	$num_topics: number of topics
	
	$dim: number of dimensions (default 2)
	
	$beta: covariance for Gaussian prior of topic coordinates (default 0.1*$num_docs)
	
	$gamma: covariance for Gaussian prior of document coordinates (default 0.1*$num_topics)
	
	$EM_iter: number of iterations for EM (default 100)
	
	$Quasi_iter: maximum iterations of Quasi-Newton (default 10)
	
	$word_vectors_file: file contains word vectors
	
	$D: dimension of word vector
	
	$rho: hyper-parameter \rho
	
	$rho_0: hyper-parameter \rho_0
	
	$data: input data
	
	$output_file: output file
	

Details:

+ This implementation needs Algorithm::LBFGS library for quasi-Newton method L-BFGS.
  The library can be downloaded at http://search.cpan.org/~laye/Algorithm-LBFGS-0.16/lib/Algorithm/LBFGS.pm.
  To install,
	
	  cpan Algorithm::LBFGS
	
+ Example of input data with 3 documents (numbers are ids of words):
	0 1 1 2 2 3 4 4 5 6 7 7 8 8 8 8 9 10 11 12 13 13 14 14 15 15 15 16
	17 18 19 20 20 21 22 23 24 25 25 25 25 25 25 26 27 27 28 29
	30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 50 51 52 53 54 54 55 56 57 58
			
+ Each line of $word_vectors_file is a word vector whose fields are separated by '\t'. The line numbers corresponding to the ids of words.

HOW TO CITE
If you use GaussianSV for your research, please cite:

	@inproceedings{gaussiansv,
	    title={Semantic Visualization for Short Texts with Word Embeddings},
	    author={Le, Tuan MV and Lauw, Hady W},
	    booktitle={International Joint Conference on Artificial Intelligence},
	    year={2017}
	}
		
	The paper can be downloaded from: http://www.hadylauw.com/publications/ijcai17a.pdf
