#!/usr/bin/env perl

# Fixes a fasta file whose deflines do not have proper
# newlines around it. Run the script with -h for help.
# Author: Lee Katz <lkatz@cdc.gov>

use warnings;
use strict;
use Getopt::Long;
use File::Basename qw/basename/;

exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help));
  die usage() if($$settings{help});

  my $lineCount=0;
  my $contigCount=0;
  while(<>){
    $lineCount++;
    s/^\s+|\s+$//g; # super-chomp, gets rid of \r also

    # Find any contigs and give the deflines a newline
    # before and after. Give a new contig name.
    s/>/\n>/g;

    # Find any places where a digit or other weird char
    # mashes up against a nucleotide.
    s/([\d_-])([ATCGN])/$1\n$2/gi;

    # Look at individual lines now that newlines have
    # been reinserted.
    for my $line(split(/\n/,$_)){ 
      $line=~s/^\s+|\s+$//g;

      # Redefine the contig name
      if($line=~/>/){
        $line=sprintf(">contig%06d",++$contigCount);
      }
      # Remove lines that are obviously fragments 
      # from a chopped defline. Those that have
      # something other than ATCGN.
      else {
        next if($line=~/[^ATCGN]/i);
      }

      # Must not be blank
      next if($line=~/^\s*$/);


      # Finally!  Print.
      print $line."\n";
    } 
  }

  # If no lines have been given to this script,
  # then print the usage.
  if($lineCount==0){
    die usage();
  }

  return 0;
}

sub usage{
  local $0=basename($0);
  "$0: fixes newlines around deflines in a fasta file.
Reads from STDIN and writes to STDOUT.

  EXAMPLE USAGE
       $0 < mashedDeflines.fasta > corrected.fasta
       tr --delete '\\n' < file.fasta | $0 > corrected.fasta
       echo '>contig1AATTGG' | $0 > corrected.fasta
  "
}
