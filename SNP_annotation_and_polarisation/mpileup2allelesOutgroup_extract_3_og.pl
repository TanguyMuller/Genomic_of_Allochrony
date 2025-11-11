#!/usr/bin/perl
use strict;
use warnings;

# Extract alleles from an mpileup file for outgroups.
# Edit $numCol_lib_min and $numCol_lib_max as needed.

my $file_mpileup = shift @ARGV; 
my $output = shift @ARGV; 

# Need adapt depends of the number of outgroups
my $numCol_lib_min = 3;  # first mpileup column to process
my $numCol_lib_max = 9;  # last mpileup column to process

extract_data_from_mpileup($file_mpileup, $output);

sub extract_data_from_mpileup {
    my ($mpileup, $output) = @_;
    my ($nb_SNP, $nb_SNPbi) = (0, 0);

    open(my $OUT, ">", $output) or die "Cannot write to $output: $!\n";

    # Need adapt with your names
    print $OUT "SNPbi\tREF\tAF_c\tAF_allele\tAF_STATUT\tGU_c\tGU_allele\tGU_STATUT\tXI_c\tXI_allele\tXI_STATUT\tNbOutgr\n";

    open(my $MPI, "<", $mpileup) or die "Cannot open $mpileup: $!\n"; 
    while (my $line = <$MPI>) {
        next if $line =~ /^#/;
        chomp $line;
        my @a_line = split(/\t/, $line);
        my $nb_outgr = 0;
        $nb_SNP++;

        my $snp = $a_line[0] . "__" . $a_line[1];
        my $ref = uc($a_line[2]);
        $ref =~ s/N/N/;
        $snp .= "\t$ref\t";

        for (my $i = $numCol_lib_min; $i <= $numCol_lib_max; $i += 3) {
            my $cov = $a_line[$i];
            my $listalleles = $a_line[$i + 1];
            $listalleles =~ s/\$//g;
            $listalleles =~ s/\^.//g;

            if ($cov == 0) {
                $snp .= "$cov\tNA\tCOVnul\t";
            } else {
                if ($cov == length($listalleles) && $listalleles !~ /[+-]/) {
                    my ($nb, $distinct_letter) = distinct_alleles($listalleles);
                    $distinct_letter =~ s/R/$ref/;
                    if ($nb == 1) {
                        $snp .= "$cov\t$distinct_letter\tUniq\t";
                    } else {
                        $snp .= "$cov\t$distinct_letter\tMult\t";
                    }
                    $nb_outgr++;
                } else {
                    $snp .= "$cov\t$listalleles\tOthers\t";
                }
            }
        }
        $snp .= $nb_outgr;
        print $OUT "$snp\n";
    }
    close $MPI;
    close $OUT;
    print "Nb de positions (input) = $nb_SNP\n";
}

sub distinct_alleles {
    my ($str) = @_;
    $str =~ s/[.,]/R/g;
    $str =~ tr/agctn/AGCTN/;
    $str =~ s/\*/D/g;

    my %h_alleles;
    $h_alleles{$_}++ for split //, $str;

    my $nb_allele = scalar keys %h_alleles;
    my $distinct_letter = join(";", map { "$_:$h_alleles{$_}" } keys %h_alleles);

    $distinct_letter =~ s/;$//;
    if ($nb_allele == 1) {
        $distinct_letter =~ /^(\w)/;
        $distinct_letter = $1;
    }
    return ($nb_allele, $distinct_letter);
}
