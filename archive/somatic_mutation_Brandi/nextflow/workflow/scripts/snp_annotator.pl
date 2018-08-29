#!/usr/bin/perl -w
#snp_annotator.pl

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
#use DBI;
require '/qbrc/home/bcantarel/seqprg/scripts/config.pl';
my %opt = ();
my $results = GetOptions (\%opt,'help|h','input|i=s','tfam|t=s',
			  'assembly|a=s','cancer|c');

my $vcf = $opt{input};
$opt{gver} = 'GRCh38.82' unless ($opt{gver});

my $prefix = $vcf;
$prefix =~ s/\.vcf.*$//;

my $runannot = "$java -Xmx10g -jar $snpeff -no-intergenic -lof -c $snpeff_config $opt{gver} $vcf";
$runannot .= " | $java -Xmx10g -jar $snpsift annotate $refdir/dbSnp.vcf.gz -";
$runannot .= " | $java -Xmx10g -jar $snpsift annotate $refdir/clinvar-latest.vcf.gz -";
$runannot .= " | $java -Xmx10g -jar $snpsift annotate $refdir/ExAC.vcf.gz -";
$runannot .= " | $java -Xmx10g -jar $snpsift caseControl -tfam $opt{tfam} -" if ($opt{tfam});
$runannot .= " | $java -Xmx10g -jar $snpsift annotate $refdir/cosmic.vcf.gz -" if ($opt{cancer});
$runannot .= " | $java -Xmx10g -jar $snpsift dbnsfp -v -db $refdir/dbNSFP.txt.gz -";
$runannot .= " | $java -Xmx10g -jar $snpsift gwasCat -db $refdir/gwas_catalog.tsv -";
$runannot .= " | bgzip > $prefix\.annot.vcf.gz\n";
system($runannot) unless (-e "$prefix\.annot.vcf.gz");

# open VCF, "gunzip -c $prefix\.annot.vcf.gz |" or die $!;
# while (my $line = <VCF>) {
#   chomp($line);
#   if ($line =~ m/^#CHROM/) {
#     my @header = split(/\t/,$line);
#     ($chrom, $pos,$id,$ref,$alt,$score,
#      $filter,$info,$format,@subjacc) = split(/\t/, $line);
#   }
#   next if ($line =~ m/^#/);
#   my ($chrom, $pos,$id,$ref,$alt,$score,
#       $filter,$annot,$format,@gts) = split(/\t/, $line);
#   next if ($ref =~ m/\./ || $alt =~ m/\./ || $alt=~ m/,X/);
#   next unless($chrom =~ m/chr\d+$|chrM|chrX/);
#   my %hash = ();
#   foreach $a (split(/;/,$annot)) {
#     my ($key,$val) = split(/=/,$a);
#     $hash{$key} = $val unless ($hash{$key});
#   }
#   my %evsannot;
#   if ($addannot{$chrom}{$pos}) {
#     %evsannot = %{$addannot{$chrom}{$pos}};
#     $hash{MAF} = $addannot{$chrom}{$pos}{MAF};
#     $hash{LOF} =~ s/\(|\)// if ($hash{LOF});
#     $hash{MAF}  = ',,' unless ($hash{MAF});
#   }
#   my @af = ();
#   my @afkeys = ('AMR_AF','AFR_AF','EUR_AF','SAS_AF','EAS_AF','MAF','dbNSFP_ExAC_AF','dbNSFP_ExAC_Adj_AF');
#   foreach $grp (@afkeys) {
#     next unless ($hash{$grp});
#     if ($grp eq 'MAF') {
#       push @af, split(/,/,$hash{$grp});
#     }
#     else {
#       if ($hash{$grp} =~ m/,/) {
# 	$hash{$grp} = 0.31;
#       }
#       push @af, $hash{$grp}*100;
#     }
#   }
#   @af = sort {$b <=> $a} @af;
#   next if ($af[0] && $af[0] > 25);
#   #$present = $dbh->selectall_arrayref(qq{select locus from varannot where chrom='$chrom' and pos=$pos});
#   my $locusid;
#   if (scalar(@{$present}) > 0) {
#     $snp = shift @{$present};
#     $locusid = $snp->[0];
#   }else {
#     %ins1 = (chrom=>$chrom,pos=>$pos,id=>$id,refal=>$ref,altal=>$alt,
# 	     eur_af=>$hash{EUR_AF},eas_af=>$hash{EAS_AF},sas_af=>$hash{SAS_AF},
# 	     afr_af=>$hash{AFR_AF},amr_af=>$hash{AMR_AF},
# 	     evs_ea_af=>$hash{dbNSFP_ESP6500_EA_AF},
# 	     evs_aa_af=>$hash{dbNSFP_ESP6500_AA_AF},exac_af=>$hash{dbNSFP_ExAC_AF},
# 	     exac_adj_af=>$hash{dbNSFP_ExAC_Adj_AF},clinvar_disease=>$hash{CLNDBN},
# 	     clinvar_sig=>$hash{CLNSIG},clinical_assoc_evs=>$evsannot{CA},
# 	     gwascat=>$hash{GWASCAT_TRAIT});
#     #$locusid = insert('varannot',\%ins1);
#   }
#   my %trxdone;
#   if ($hash{ANN}) {
#     foreach $trx (split(/,/, $hash{ANN})) {
#       my ($allele,$effect,$impact,$gene,$geneid,$feature,
# 	  $featureid,$biotype,$rank,$codon,$aa,$pos_dna,$len_cdna,
# 	  $cds_pos,$cds_len,$aapos,$aalen,$distance,$err) = split(/\|/,$trx);
#       my $key = join("|",$effect,$gene,$aa);
#       next if ($trxdone{$key});
#       if ($aa) {
# 	#$effpres = $dbh->selectall_arrayref(qq{select vareffid from vareff where chrom='$chrom' and pos=$pos and effect='$effect' and prot='$aa'});
#       }else {
# 	#$effpres = $dbh->selectall_arrayref(qq{select vareffid from vareff where chrom='$chrom' and pos=$pos and effect='$effect'});
#       }
#       next if (scalar(@{$effpres}) > 0);
#       next unless ($impact =~ m/HIGH|MODERATE/ || $effect =~ m/TF|regulatory/);
#       next if ($impact =~ m/MODIFIER/ && length($info{REF}) + length($info{ALT}) > 2);
#       my $exonnum = '';
#       ($exonnum, $exontot) = split(/\//,$rank) if ($rank);
#       $trxdone{$key} = 1;
#       my %ins2 = (locus=>$locusid,chrom=>$chrom,pos=>$pos,effect=>$effect,impact=>$impact,
# 		  symbol=>$gene,ensembl=>$geneid,prot=>$aa,dna=>$codon,lof=>$hash{LOF},
# 		  nmd=>$hash{NMD},gerp_rs=>$hash{dbNSFP_GERP___RS},sift=>$hash{dbNSFP_SIFT_pred},
# 		  polyphen2_hvar=>$hash{dbNSFP_Polyphen2_HVAR_pred},exon_num=>$exonnum,
# 		  phastcons100=>$hash{dbNSFP_phastCons100way_vertebrate},
# 		  mutationtaster=>$hash{dbNSFP_MutationTaster_pred},
# 		  interpro=>$hash{dbNSFP_Interpro_domain});
#       #$effid = insert('vareff',\%ins2);
#     }
#   }
#   @posalls = ($ref,split(/,/,$alt));
#   $ac = 0;
#   $an = 0;
#   @deschead = split(":",$format);
#  F1:foreach $subjid (@subjacc) {
#     my $allele_info = shift @gts;
#     @ainfo = split(/:/, $allele_info);
#     my %gtinfo = ();
#     foreach $k (0..$#deschead) {
#       $gtinfo{$deschead[$k]} = $ainfo[$k];
#     }
#     my ($x,$y) = (split(/\//, $gtinfo{GT}));
#     next F1 if ($x eq '.' || $y eq '.');
#     if ($gtinfo{AD}) {
#       @ar= split(/,/,$gtinfo{AD});
#       $gtinfo{DP} = sum(@ar) unless $gtinfo{DP};
#       $gtinfo{RO} = shift @ar;
#       $gtinfo{AO} = join(",",@ar);
#     }
#     if ($gtinfo{DV} && $gtinfo{DP}) {
# 	$gtinfo{AO} = $gtinfo{DV};
# 	$gtinfo{RO} = $gtinfo{DP} - $gtinfo{DV};
#     }elsif ($gtinfo{DPR}) {
# 	@ar= split(/,/,$gtinfo{DPR});
# 	$gtinfo{DP} = sum(@ar) unless $gtinfo{DP};
# 	$gtinfo{RO} = shift @ar;
# 	$gtinfo{AO} = join(",",@ar);
#     }
#     if ($gtinfo{VR} || $gtinfo{RR} || $gtinfo{NV} || $gtinfo{NR}) {
# 	$gtinfo{AO} = $gtinfo{VR} if ($gtinfo{VR});
# 	$gtinfo{RO} = $gtinfo{RR} if ($gtinfo{RR});
# 	$gtinfo{AO} = $gtinfo{NV} if ($gtinfo{NV});
# 	$gtinfo{RO} = $gtinfo{NR} if ($gtinfo{NR});
#     }
#     $gtinfo{AO} = 0 if ($gtinfo{DP} == $gtinfo{AO});
#     next F1 if ($gtinfo{DP} && $gtinfo{DP} < 10);
#     next F1 if ($gtinfo{GQ} && $gtinfo{GQ} < 5);
#     my %ins3 = (sample_name=>$subjid,locus=>$locusid,chrom=>$chrom,pos=>$pos,
# 		quality=>$score,filter=>$filter,gt=>$gtinfo{GT},dp=>$gtinfo{DP},
# 		alt_obs=>$gtinfo{AO},ref_obs=>$gtinfo{RO});
#     #$samsnpid = insert('sample_snp',\%ins3);
#   }
# }

