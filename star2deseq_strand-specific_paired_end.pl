#!usr/bin/perl -w
use strict;
use warnings;
use Data::Dumper;

my $file1= "sample_group_folder_details.txt";
open IN, $file1 or die;
my @sample_details = <IN>;

my %sample_database;
foreach my $line1 (@sample_details)
  {
    chomp $line1;
    my @arr1 = split("\t", $line1);
    $sample_database{$arr1[0]}{"sample"} = $arr1[1];
    $sample_database{$arr1[0]}{"group"} = $arr1[2];
  }
# print Dumper \%sample_database;

open OUT, ">counttable.txt" or die;
my @read_files = glob('*/STAR_Alignment*/ReadsPerGene.out.tab');
# print Dumper \@read_files;
my %mega_database;
foreach my $f (@read_files)
  {
    chomp $f;
    # print $f,"\n";
    my @arr2 = split("\/", $f);
    # print Dumper \@arr2;
    my $folder = $arr2[6];
    chomp $folder;
    my $sample = $sample_database{$folder}{"sample"};
    my $group = $sample_database{$folder}{"group"};
    #######Open file and process########
    open IN1, $f or die;
    my @file_content = <IN1>;
    foreach my $line2 (@file_content)
      {
        chomp $line2;
        next if $line2 =~ /^(N_)/;
        # print $line2, "\n";
        my @arr3 = split("\t", $line2);
        my $gene = $arr3[0];
        my $count = $arr3[3];
        $mega_database{$gene}{$sample} = $count;

      }
  }
# print OUT Dumper \%mega_database;
my @final_output;
my @total_gene_set;
foreach my $gene (sort keys %mega_database)
  {
    push(@total_gene_set, $gene);
  }
# print OUT Dumper \@total_gene_set;
push (@final_output, "genes", "\t");
my @total_sample;
foreach my $sample (sort keys %{$mega_database{$total_gene_set[0]}})
    {
      push(@total_sample, $sample);
      push(@final_output, $sample, "\t");
    }
push (@final_output, "\n");
# print OUT Dumper @final_output;
# foreach my $a (@final_output)
#   {
#     print OUT $a;
#   }
foreach my $gene (sort keys %mega_database)
  {
    push (@final_output, $gene, "\t");
    foreach my $sample (sort keys %{$mega_database{$gene}})
          {
            my $count = $mega_database{$gene}{$sample};
            push(@final_output, $count,"\t");
          }
    push(@final_output, "\n");
  }
foreach my $a (@final_output)
  {
    # chomp $a;
    print OUT $a;
  }

####Coldata###
open OUT2, ">coldata.txt" or die;
print OUT2 "\tcondition\ttype\n";
foreach my $line3 (@sample_details)
  {
    chomp $line3;
    my @arr4 = split("\t", $line3);
    print OUT2 "$arr4[1]\t$arr4[2]\tpaired-end\n";
  }
