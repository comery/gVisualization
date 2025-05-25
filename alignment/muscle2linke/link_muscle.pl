#!/usr/bin/perl -w 
###################################################
#this perl is to change the format of muscle output
#both DNA and PRO seq can be accepted!
###################################################
=head1 name

	link_muscle.pl
	
=head1 Author
	
	Chentao YANG 
	yangchentao@genomics.cn

=head1 Description
	
	This program is to change the -clw format runed 
	by muscle, and provide a easy way to get statistics
	of match & mismatch information. In addition, you
	can also use -trim/-primer to get more infos.

=head1 Version
	
	Version: 1.0,	Date:2015-1-12
	Version: 2.0,	Date:2016-10-27
	can compare hiseq sequence with sanger sequence.
	Version: 3.0,	Date:2018-06-27
	input format add 'fasta', and fix bug.

=head1 Usage
	
	perl $0  *.clw|*.fasta [option]
	--primer <str> primer string
	--trim whether you want to do triming
	--del <number>	how many nucleic acids you want to trim
	--tag <str> which sequence you want to set as standard to do triming,here you just give id [>id].
	--sanger give some alignment information of degenerate base, but it just works when sequences total number is exactly 2!
	--out <str> output file including special symbles.
	
=head1 Exmple

	perl link_muscle.pl 2.clw -o 2.link
	perl link_muscle.pl 2.clw -trim -del 30 -sanger -out 2.link

	
=cut


use strict;
use List::Util qw(max);
use Getopt::Long;
my ($trim,$del,$tag,$primer,$sanger,$output,$help);
GetOptions(
	"trim!"=>\$trim,
    "del:i"=>\$del,
    "tag=s"=>\$tag,
    "primer!"=>\$primer,
    "sanger!"=>\$sanger,
    "out|o:s"=>\$output,
    "help!"=>\$help
    );

#print "$trim,$del,$tag,$primer,$sanger,$output\n";
#die "Usage\n\tperl $0 input.clw trimed_number \n" unless (@ARGV == 2);
die `pod2text $0` if ($help || @ARGV==0);
$del ||= 50;
$output ||= "$ARGV[0].link";

open IN,"$ARGV[0]" or die "$!!";
open OUT,">$output";

my ($xingxing,%hash,$id,%id_len,$seq,%code);

%code= (
			'-' => 0,
			'.' => 0,
			'A' => 2,
			'T' => 3,
			'C' => 5,
			'G' => 7,
			'R' => 14,		# A|G
			'Y' => 15,		# C|T
			'M' => 20,		# A|C
			'K' => 21,		# G|T
			'S' => 35,		# G|C
			'W' => 24,		# A|T here is not 6 for setting a cufoff to be easy to know jianbing base
			'H' => 30,		# A|T|C
			'B' => 105,		# G|T|C
			'V' => 70,		# G|A|C
			'D' => 42,		# G|A|T
			'N' => 210		# A|G|C|T

);

my $check_file_type = `head -1 $ARGV[0]`;
my $format;
if ($check_file_type =~ /^>/){
	print OUT "#input file format: FASTA\n";
	$format = 1;
}elsif ($check_file_type =~ /^MUSCLE/){
	print OUT "#input file format: CLUSTALW\n";
	$format = 0;
}else{
	print STDERR "can not regonize input file!\n";
	exit;
}

my $max_id_len = 0;
if ($format){
	$/=">";<IN>;$/="\n";
	while ($id =<IN>){
		chomp $id;
		$/=">";
		chomp(my $seq = <IN>);
		$seq =~ s/\n//g;
		$seq = uc($seq);
		$/="\n";
		my $ll = length $id;
		$id_len{$id} = $ll;
		$max_id_len = $ll if ($ll > $max_id_len);
		if(exists $hash{$id}){
			$hash{$id}="$hash{$id}"."$seq";
		}else{
			$hash{$id}=$seq;
		}
	}

}else{
	while (<IN>){
		$_=~s/^\s+$//g;
		chomp;
		if ($_ eq ""||$_=~/^MUSCLE/|| $_=~/\*/) {
			next;
		}else{
			$id=(split /\s+/,$_)[0];
			$seq=(split /\s+/,$_)[-1];
			my $ll = length $id;
			$id_len{$id} = $ll;
			$max_id_len = $ll if ($ll > $max_id_len);
			if(exists $hash{$id}){
				$hash{$id}="$hash{$id}"."$seq";
			}else{
				$hash{$id}=$seq;
			}
		}
	}
}

my ($key,$value,@seq,$blank);

while (($key,$value)=each %hash){
	push @seq,$value;
	$blank = " " x ($max_id_len - $id_len{$key} + 8 );
	print OUT "$key"."$blank"."$value\n";
}
print OUT ">>consensus"; # this word's length is 11;
print OUT " " x ($max_id_len + 8 - 11);

my $length=length($seq[0]);
my ($y,@array,$symble);

for($y=0;$y<$length;$y++){
	foreach my $align(@seq){
		my $aa=substr $align,$y,1;
		push @array,$aa;
	}
	my $tmp = &judge(@array);
	$symble .= "$tmp";
	@array=();
}
print OUT "$symble\n";

#-------BASIC MODEL----------------------
my $ll = length ($symble);

my $mismatch = ($symble =~ tr/\|/\|/);
my $conserved = ($symble =~ tr/\*/\*/);
my $gap_whole = ($symble =~ tr/-/-/);
my $rate = $conserved/$ll;
print OUT "\n----BASIC INFO----\n";
print OUT "whole_alignment_length\t\t$ll\nConserved\t$conserved\n";
print OUT "gap_whole\t$gap_whole\n";
print OUT "mismatch_whole\t$mismatch\n";
print OUT "conserved_rate_whole:";
printf OUT "%.3f\n",$rate;


#------TRIME MODEL---------------------
my ($care,$gap);
if ($trim) {
	if ($tag) {
		my $standard;
		foreach my $j(keys %hash){
			$standard = $j if ($j =~ /$tag/);
		}
		my ($sub_l,$sub_r) = &trim($hash{$standard});
		my $trimed_ll = $ll - $sub_l - $sub_r;
		my $end_care = $sub_l + $trimed_ll;
		$care = substr($symble,$sub_l,$trimed_ll);
		$gap = ($care =~ tr/-/-/);
		my $same = ($care =~ s/\*/\*/g);
		my $mis = ($care =~ tr/\|/\|/);

		print OUT "\n----TRIME INFO----\n";
		print OUT "The information of region [$sub_l-$end_care] is showed:\n";
		print OUT "care_length\t$trimed_ll\n";
		print OUT "conserved_care\t$same\n";
		print OUT "mismatch_care\t$mis\n";
		print OUT "gap_care\t$gap\n";
		my $ratio = $same/$trimed_ll;
		print OUT "conserved_rate_care\t";
		printf OUT "%.3f\n",$ratio;
	}else {
		print STDOUT "skipping because of no tag! ";
		last;
	}
	
}

close IN;

sub trim {
	my $str = shift;
	my $base_s = 0;
	my $count_s = 0;
	while ($base_s < $del) {
		my $tmp = substr($str,$count_s,1);
		$base_s ++ if ($tmp =~ /[A|T|C|G]/);
		$count_s ++;
	}

	my $rev=reverse $str;
	my $base_e = 0;
	my $count_e = 0;
	while ($base_e < $del) {
		my $tmp = substr($rev,$count_e,1);
		$base_e ++ if ($tmp =~ /[A|T|C|G]/);
		$count_e ++;
	}

	return ($count_s,$count_e);
}

#--------PRIMER MODEL-------------------
my $case_p;
if ($primer) {
	if ($trim && $care) {
		$case_p = $care;
	}else{
		$case_p = $symble;
	}
	my $repeat=&max_repeat($case_p);
	my $blocks = &primer ($case_p);
#	my $xing_lens = &index($case_p);
	print OUT "\n----PRIMER INFO----\n";
	print OUT "max_repeat\t$repeat\n";
	print OUT "potential_primer_blocks\t$blocks\n";
#	print OUT "xing_lens\t$xing_lens\n";
}else {
	print STDOUT "Skipping PRIMER MODEL...\n";
}
	
sub max_repeat{
my $lens=shift;
my $long=length ($lens);
my  $j=0;my @count;my $i;
for ($i=0;$i<=$long-1;$i++){
	my $base=substr ($lens,$i,1);
	if (($base eq "*") ||( $base eq ".")){
		$j++;
	}else{
		push @count,$j;
		$j=0;
	}
}
my $max_repeat=max(@count);
return $max_repeat;
}	

sub primer {
my $str = shift;
my $len = length($str);
my $count;
my $blocks = 0;
my $n;
for ($n = 0;$n<=$len-21;$n+=20){
	my $temp_str = substr $str,$n,20;
	$count = ($temp_str=~tr/*/*/);
	if ($count >=12){
		$blocks++;
	}
}
return $blocks;
}


#------SANGER MODEL---------------------
my $case_s;
if ($sanger && keys %hash == 2 ) {
	if ($trim && $care) {
		$case_s = $care;
	}else{
		$case_s = $symble;
	}
		my $case_ll = length $case_s;
		my $same = ($case_s =~ s/\*/\*/g);
		my $jb_care_right = ($case_s =~ tr/\./\./);
		my $jb_care_wrong = ($case_s =~ tr/x/x/) ;
		my $gap_care = ($case_s =~ tr/-/-/) ;

		print OUT "\n----SANGER INFO----\n";
		print OUT "jb_care_right\t$jb_care_right\n";
		print OUT "jb_care_wrong\t$jb_care_wrong\n";
		print OUT "gap_care\t$gap_care\n";
		my $ratio = ($same+$jb_care_right)/$case_ll;
		print OUT "conserved_rate_real\t";
		printf OUT "%.3f\n",$ratio;
		print OUT "\n# print some usful info:";
		print OUT "there are five alignment pattern:\n";
		print OUT "\"-\" means gap\n";
		print OUT "\"*\" means same nucleic acid\n";
		print OUT "\"|\" means mismatch, any of [A|T|C|G]\n";
		print OUT "\".\" means one is [R|Y|M|K|S|W|H|B|V|D|N] And another is which before one should contain\n";
		print OUT "\"x\" means one is [R|Y|M|K|S|W|H|B|V|D|N] But another is which before one dosen't contain\n";
}elsif ($sanger && keys %hash != 2) {
	die "if you want use -sanger mode, must be sure input sequences number = 2!\n";
}else{
	print STDOUT "Skipping SANGER MODEL...\n";
}

sub judge {
	my @aa = @_;
	if (@aa == 2) {
		my $a = $aa[0];
		my $b = $aa[1];

		if ($code{$a} == $code{$b} && $code{$a}*$code{$b} != 0) {
			return "*";
		}elsif ($code{$a}*$code{$b} == 0 ){
			return "-";
		}elsif ($code{$a} != $code{$b} && ($code{$a} <= 7 && $code{$b} <= 7)){
			return "|";
		}elsif ($code{$a} + $code{$b} >14 && &cont($a,$b)) {
			return ".";
		}elsif ($code{$a} + $code{$b} >14 && !&cont($a,$b)){
			return "x";
		}else {
			return "U";
		}
	}else {
		my %count;
		my @uniq_aa = grep { ++$count{ $_ } < 2; } @array;
		my $mulit_i = 1;
		foreach my $i (@aa){
			$mulit_i = $code{$i} * $mulit_i;
		}
		if ($mulit_i == 0 ){
			return "-";
		}elsif(@uniq_aa == 1 && $uniq_aa[0] =~ /[A|T|C|G]/){
			return "*";
		}else {
			return "|";
		}
	}

}

 
sub cont {
	my @aa = @_;
	my $a = $aa[0];
	my $b = $aa[1];
	if ($code{$a}/$code{$b} >1 && $code{$a} % $code{$b} ==0) {
		return 1;
	}elsif ($code{$b}/$code{$a} >1 && $code{$b} % $code{$a} ==0){
		return 1;
	}else {
		return 0;
	}
}
