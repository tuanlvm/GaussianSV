# Copyright 2017 Singapore Management University (SMU). All Rights Reserved. 
#
# Permission to use, copy, modify and distribute this software and 
# its documentation for purposes of research, teaching and general
# academic pursuits, without fee and without a signed licensing
# agreement, is hereby granted, provided that the above copyright
# statement, this paragraph and the following paragraph on disclaimer
# appear in all copies, modifications, and distributions.  Contact
# Singapore Management University, Intellectual Property Management
# Office at iie@smu.edu.sg, for commercial licensing opportunities.
#
# This software is provided by the copyright holder and creator "as is"
# and any express or implied warranties, including, but not Limited to,
# the implied warranties of merchantability and fitness for a particular 
# purpose are disclaimed.  In no event shall SMU or the creator be 
# liable for any direct, indirect, incidental, special, exemplary or 
# consequential damages, however caused arising in any way out of the
# use of this software.


use warnings;
use strict;
use Data::Dumper;
use Math::Trig;
use Algorithm::LBFGS;

use Getopt::Long;

our $num_topics = 10;
our $data;
our $dim = 2;
our $num_iter = 100;
our $num_quasi_iter = 10;
our $output_file;
our $beta;
our $gamma;
our $word_vectors_file;
our $D=300; #The dim of the vectors
our $rho_z_const = 100;
our $rho_0 = 10000;

GetOptions ("num_topics=i" => \$num_topics,
			"dim=i" => \$dim,
            "beta=f"     => \$beta,
            "gamma=f"   => \$gamma,
            "EM_iter=i" => \$num_iter,
            "Quasi_iter=i" => \$num_quasi_iter,
            "data=s" => \$data,
            "word_vectors=s" => \$word_vectors_file,
            "word_vector_dim=i" => \$D,
            "rho=f" => \$rho_z_const,
            "rho_0=f" => \$rho_0,
            "output_file=s" => \$output_file)
  or die("Error in command line arguments\n");

open(my $output, '>', $output_file) or die "Could not open file '$output_file' $!";

$gamma  = $gamma ? $gamma : 0.1 * $num_topics;

our $num_docs = 0;
our @input_lines;
my %Vocabulary;
open( F, "$data" ) or die "Couldn't load the data $data $!";
while (<F>) {
	chomp;
	my $line = $_;
	push(@input_lines, $line);
	$num_docs++;
}
close(F);
$beta = $beta ? $beta : 0.1 * $num_docs;

print $output "num_topics: ";
print $output $num_topics;
print $output "\n";

print $output "dim: ";
print $output $dim;
print $output "\n";

print $output "beta: ";
print $output $beta;
print $output "\n";

print $output "gamma: ";
print $output $gamma;
print $output "\n";

print $output "EM_iter: ";
print $output $num_iter;
print $output "\n";

print $output "Quasi_iter: ";
print $output $num_quasi_iter;
print $output "\n";

print $output "word vectors file: ";
print $output $word_vectors_file;
print $output "\n";

print $output "rho_z_const: ";
print $output $rho_z_const;
print $output "\n";

print $output "rho_0: ";
print $output $rho_0;
print $output "\n";

our @word_vectors;#row-based
#load word vectors
open( F,
	"$word_vectors_file"
	) or die "Couldn't open $word_vectors_file $!";
my $num_word = 0;
while (<F>) {
	chomp;
	my $line = $_;
	my @tokens = split( /\s/, $line );
	my @position;
	for ( my $d = 0 ; $d < $D ; $d++ ) {
		$position[$d]  = $tokens[$d];
	}
	$word_vectors[$num_word] = \@position;
	$num_word++;
}
close(F);

our @mu_0_list;
for(my $d=0;$d<$D;$d++){
	my $sum = 0;
	for(my $k = 0; $k < $num_word; $k++ ){
		$sum += $word_vectors[$k][$d];
	}
	$mu_0_list[$d] = $sum/$num_word;
}

our %gaussian_zw_list;
our %mu_z_list;
#our %inverse_sigma_z_list;

for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
	my %list;
	my @array;
	my %lst;
	for(my $j=0;$j<$num_word;$j++){
		$array[$j] = 0;
		my @a;
		my @b;
		for(my $d=0;$d<$D;$d++){
			$a[$d] = 0;
			$b[$d] = 0;
		}
		$list{$j} = \@a;
		$lst{$j} = \@b;
	}
	$gaussian_zw_list{$i} = \@array;

	my @arr_inverse;
	my @arr_mu;
	for(my $j=0;$j<$D;$j++){
		$arr_mu[$j] = &gaussian_rand * 0.01 + $mu_0_list[$j];
	}
	$mu_z_list{$i} = \@arr_mu;
}
	
my %all_tokens;
my %all_tokens_count;
my @doc_length;
my $c = 0;
foreach (@input_lines) {
	my $line = $_;
	my @tokens = split( /\s/, $line );
	$doc_length[$c] = @tokens;
	my @unique_tokens;
	my @count_unique_tokens;
	my $count = 0;
	my %hash_w;
	for ( my $i = 0 ; $i < @tokens ; $i++ ) {
		if(!defined $Vocabulary{$tokens[$i]}) {
			$Vocabulary{$tokens[$i]} = 1;
		}
		if(defined $hash_w{$tokens[$i]}) {
			$hash_w{$tokens[$i]} = $hash_w{$tokens[$i]} + 1;
		}
		else{
			$hash_w{$tokens[$i]} = 1;
		}
	}
	foreach my $name (keys %hash_w) {
   		$unique_tokens[$count] = $name;
   		$count_unique_tokens[$count] = $hash_w{$name};
   		$count++;
	}
	$all_tokens{$c} = \@unique_tokens;
	$all_tokens_count{$c} = \@count_unique_tokens;
	$c++;
}

our $wordsize = (keys %Vocabulary);

# init parameters
our %phi;
our %xai;

our %thetaupdate_numerator;


# init phi
print $output "init phi \n";
for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
	my @position;
	for ( my $d = 0 ; $d < $dim ; $d++ ) {
		$position[$d]  = &gaussian_rand * 0.01;
		print $output $position[$d];
		print $output "\n";
	}
	print $output "---------------\n";
	$phi{$i}     = \@position;
}

# init xai
print $output "init xai \n";

for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
	my @position;
	for ( my $d = 0 ; $d < $dim ; $d++ ) {
		$position[$d]  = &gaussian_rand * 0.01;
		print $output $position[$d];
		print $output "\n";
	}
	print $output "---------------\n";
	$xai{$i}     = \@position;
}

our %pr_zgx;
our %pr_zgnm;
our %sum_zgnm_overw;

my @phi_gradient;
my @xai_gradient;

my $quasi   = Algorithm::LBFGS->new;
$quasi->set_param(max_iterations => $num_quasi_iter);
my $eval_short = sub {
	my $x = shift;
	my $sum   = 0;
	my $docid = 0;
	my %grad;
	my @gradient;
	my @numerators;
	for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
		my @arr;
		for ( my $k = 0 ; $k < $dim ; $k++ ) {
			my $index = $j * $dim + $k;
			$arr[$k] = - $beta * $x->[$index];
			$gradient[$index] = 0;
		}
		$grad{$j} = \@arr;
	}
	for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
		my @probs;
		my $denominator = 0;
		my $offset_doc = $num_topics * $dim + $i*$dim;
		for ( my $k = 0 ; $k < $num_topics ; $k++ ) {
			my $d = 0;
			my $offset_topic = $k * $dim;
			for ( my $l = 0 ; $l < $dim ; $l++ ) {
				$d += ($x->[$offset_doc + $l] - $x->[$offset_topic + $l])**2;
			}
			$numerators[$k] = exp((-0.5) * $d );
			$denominator += $numerators[$k];
		}
		for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
			$probs[$j] = $numerators[$j] / $denominator;
		}
		$pr_zgx{$i} = \@probs;
		
		#compute Q
		$docid = $i;
		my @tokens = values($all_tokens{$docid});
		my @count_tokens = values($all_tokens_count{$docid});
		 
		my $p_zgnm = $pr_zgnm{$docid};
		my $p_zpx  = $pr_zgx{$docid};
		for ( my $m = 0 ; $m < @tokens ; $m++ ) {
			my $p_znm = $p_zgnm->{$m};
			my $t = 0;
			for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
				$t += 
				  $p_znm->[$j] *
				  (log( $p_zpx->[$j]));
			}
			$sum += $count_tokens[$m] * $t;
		}
		
		my $length = 0;
		for ( my $k = 0 ; $k < $dim ; $k++ ) {
			$length += $x->[$offset_doc + $k]**2;
		}
		$sum += (-$gamma/2) * $length;
		
		my @gra;
		for ( my $k = 0 ; $k < $dim ; $k++ ) {
			$gra[$k] = 0;
		}
		#gradient wrt phi
		for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
			my $offset_topic = $j * $dim;
			my $p_zgnm = $pr_zgnm{$docid};

			my $euclid = 0;
			for ( my $l = 0 ; $l < $dim ; $l++ ) {
				$euclid += ($x->[$offset_doc + $l] - $x->[$offset_topic + $l])**2;
			}			

			my $p_zpx = $pr_zgx{$docid}->[$j];
			for ( my $k = 0 ; $k < $dim ; $k++ ) {
				my $t;
				$t = ( $doc_length[$docid] * $p_zpx - $sum_zgnm_overw{$docid}->[$j] ) * ( $x->[$offset_topic + $k] - $x->[$offset_doc + $k] );
				$grad{$j}->[$k] +=  $t;
				#wrt xai
				$gra[$k] += -$t;
			}
		}
		for ( my $k = 0 ; $k < $dim ; $k++ ) {
			my $diff = $gra[$k] - $gamma * $x->[$offset_doc + $k];
			$gradient[$offset_doc + $k] = -$diff;
		}
	}
	for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
		$sum += ( -$beta / 2 ) * ( $x->[$j*$dim + 0]**2 + $x->[$j*$dim + 1]**2);
	}
	
	my $f = -$sum;
	print ($f);
	print ("\n");
	
	print $output $f;
	print $output "\n";
	
	for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
		for ( my $k = 0 ; $k < $dim ; $k++ ) {
			$gradient[$j * $dim + $k] = -$grad{$j}->[$k];
		}
	}
	return ( $f, \@gradient );
};

#EM
for ( my $i = 0 ; $i < $num_iter ; $i++ ) {
	print($i);
	print("\n");
	
	print $output $i;
	print $output "\n";
	
	&estep();
	&mstep();
}


# output to file
print $output "\n-------Xai-------\n";
for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
	for ( my $d = 0 ; $d < $dim ; $d++ ) {
		if($d==0){
			print $output $xai{$i}[$d];
		}
		else{
			print $output "\t";
			print $output $xai{$i}[$d];
		}
	}
	print $output "\n";
}

print $output "\n-------Phi-------\n";
for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
	for ( my $d = 0 ; $d < $dim ; $d++ ) {
		if($d==0){
			print $output $phi{$i}[$d];
		}
		else{
			print $output "\t";
			print $output $phi{$i}[$d];
		}
	}
}

print $output "\n-------Verbose-------\n";
my $docid = 0;
print $output "zgivend \n";
foreach (@input_lines) {
	my $p_zgx = $pr_zgx{$docid};
	for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
		print $output $p_zgx->[$j];
		print $output "\n";
	}
	print $output "-------------\n";
	$docid++;
}

print $output "Phi\n";
print $output Dumper(%phi);
print $output "Xai\n";
print $output Dumper(%xai);
print $output "Mu_z\n";
print $output Dumper(%mu_z_list);

sub mstep {
	#compute mu_z
	for(my $t=0;$t<$num_topics;$t++){
		for(my $d=0;$d<$D;$d++){
			my $denominator = 0;
			my $numerator = 0;
			for(my $i=0;$i<$num_docs;$i++){
				my @tokens = values($all_tokens{$i});
				my @count_tokens = values($all_tokens_count{$i});
		 
				my $p_zgnm = $pr_zgnm{$i};
				for ( my $m = 0 ; $m < @tokens ; $m++ ) {
					my $p_znm = $p_zgnm->{$m};
					my $tmp = $p_znm->[$t] * $count_tokens[$m] * $rho_z_const;
					$numerator += $tmp * $word_vectors[$tokens[$m]][$d];
					$denominator += $tmp;
				}
			}
			$numerator += $rho_0 * $mu_0_list[$d];
			$denominator += $rho_0;
			$mu_z_list{$t}->[$d] = $numerator/$denominator;
		}
	}
	
	my @x0;
	for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
		for ( my $d = 0 ; $d < $dim ; $d++ ) {
			$x0[$i * $dim + $d] = $phi{$i}[$d];
		}
	}
	for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
		for ( my $d = 0 ; $d < $dim ; $d++ ) {
			$x0[$num_topics * $dim + $i * $dim + $d] = $xai{$i}[$d];
		}
	}
	
	print("start quasi\n");
	my $x;
	$x = $quasi->fmin($eval_short, \@x0);
	print("end quasi\n");
	for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
		for ( my $d = 0 ; $d < $dim ; $d++ ) {
			$phi{$i}[$d] = $x->[$i * $dim + $d];
		}
	}
	for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
		for ( my $d = 0 ; $d < $dim ; $d++ ) {
			$xai{$i}[$d] = $x->[$num_topics * $dim + $i * $dim + $d];
		}
	}
	
	return;
}

sub estep {
	my %distances;
	for ( my $i = 0 ; $i < $num_topics ; $i++ ) {
		for ( my $j = 0 ; $j < $num_word ; $j++ ) {
			my $sum_tmp = 0;
			for(my $d=0;$d<$D;$d++){
				my $tmp = $word_vectors[$j][$d] -  $mu_z_list{$i}->[$d];
				$sum_tmp += ($tmp**2)*$rho_z_const;
			}
			my $t = exp(-0.5*$sum_tmp);
			$gaussian_zw_list{$i}->[$j] = $t;
		}
	}
	
	for ( my $i = 0 ; $i < $num_docs ; $i++ ) {
		my @probs;
		my $denominator = 0;
		my @numerators;
		for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
			my $euclid = 0;
			for ( my $k = 0 ; $k < $dim ; $k++ ) {
				$euclid += ($xai{$i}->[$k]-$phi{$j}->[$k])**2;
			}
			$numerators[$j] = exp( (-0.5) * $euclid );
			$denominator += $numerators[$j];
		}
		for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
			my $prob = $numerators[$j]/$denominator;
			$probs[$j] = $prob;
		}
		$pr_zgx{$i} = \@probs;
	}

	for ( my $docid = 0 ; $docid < $num_docs ; $docid++ ) {
		my @tokens = values($all_tokens{$docid});
		my @count_tokens = values($all_tokens_count{$docid});
		
		my %probs_znm;
		my @sum_zpnm;
		for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
			$sum_zpnm[$j] = 0;
		}
		my $p_zpx = $pr_zgx{$docid};
		for ( my $i = 0 ; $i < @tokens ; $i++ ) {
			my @probs;
			my $denominator = 0;
			my @numerator;
			for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
				my $t = $p_zpx->[$j] * $gaussian_zw_list{$j}->[$tokens[$i]];
				$denominator += $t;
				$numerator[$j] = $t;
			}
			for ( my $j = 0 ; $j < $num_topics ; $j++ ) {
				my $p = $numerator[$j]/$denominator;
				$probs[$j] = $p;
				my $temp = $count_tokens[$i]*$p;
				$sum_zpnm[$j] += $temp;
			}
			$probs_znm{$i} = \@probs;
		}
		$pr_zgnm{$docid} = \%probs_znm;
		$sum_zgnm_overw{$docid} = \@sum_zpnm;
	}
	return;
}

sub gaussian_rand {
	#This function is from Perl Cookbook, 2nd Edition (http://www.oreilly.com/catalog/perlckbk2/)
	
	my ( $u1, $u2 );    # uniformly distributed random numbers
	my $w;              # variance, then a weight
	my ( $g1, $g2 );    # gaussian-distributed numbers

	do {
		$u1 = 2 * rand() - 1;
		$u2 = 2 * rand() - 1;
		$w  = $u1 * $u1 + $u2 * $u2;
	} while ( $w >= 1 );

	$w  = sqrt( ( -2 * log($w) ) / $w );
	$g2 = $u1 * $w;
	$g1 = $u2 * $w;

	# return both if wanted, else just one
	return wantarray ? ( $g1, $g2 ) : $g1;
};
