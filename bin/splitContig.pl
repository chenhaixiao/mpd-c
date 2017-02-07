#!/usr/bin/env perl
# Name:           splitContig.pl
# Date Created:   2017-02-06T11:15:23
# Date Modified:  2017-02-06T11:15:23
# By:             TS Wingo
#
# Description:    Takes a single fasta file and splits it into multiple
#                 individual fasta files defined by the header, e.g. '>entry1',
#                 and it creates an "in" file for the index_genome program.
#                 Use the "in" file like so `index_genome < PanTig.in`
#

use 5.10.0;
use warnings;
use strict;

# common
use Getopt::Long;
use IO::Zlib;
use Path::Tiny;

# variables
my ( $fileName, $outExt );

# get options
die "Usage: $0 -f <fileName> -o <outExt>\n"
  unless GetOptions(
  'file|f=s' => \$fileName,
  'out|o=s'  => \$outExt,
  )
  and $fileName
  and $outExt;

###############################################################################
# main
###############################################################################

my ( $namesAref, $sizeOfContigHref ) = Split_fasta($fileName);
WriteInFile( $outExt, $namesAref, $sizeOfContigHref );

###############################################################################
# subroutines
###############################################################################

sub WriteInFile {
  my ( $out_ext, $namesAref, $sizeOfContigHref ) = @_;

  my @sSizes = sort { $b <=> $a } map { $sizeOfContigHref->{$_} } @$namesAref;
  my $maxSize = $sSizes[0];

  my $fh = path("$out_ext.in")->filehandle(">");
  say {$fh} join "\n", ( "D", "in.index_genome.in", scalar @$namesAref );
  say {$fh} join "\n", @$namesAref;
  say {$fh} $out_ext;
  say {$fh} ( $maxSize + 1 );
  say {$fh} "14";
}

sub Split_fasta {
  my $file = shift;

  my @names;

  my $fh = mustOpen($file);
  my ( $outFh, $fileName );
  my %sizeOfContig;

  while (<$fh>) {
    chomp $_;
    if ( $_ =~ m/\>([\w\d]+)/ ) {
      my $name = $1;
      $fileName = "$name.fa";
      say $name;
      $outFh = path($fileName)->filehandle(">");
      push @names, $fileName;
      say {$outFh} ">$name";
      next;
    }
    $sizeOfContig{$fileName} += length $_;
    say {$outFh} $_;
  }
  if ( !@names ) {
    say "Did not encounter any fasta entries, which should start like '>entry1'\n";
    exit(1);
  }
  return ( \@names, \%sizeOfContig );
}

sub mustOpen {
  my $file = shift;

  if ( $file =~ m/\.gz\z/xm ) {
    my $fh = new IO::Zlib;
    if ( !-e $file ) {
      die "Error: Cannot open '$file' - file does not exist.\n";
    }
    elsif ( !$fh->open( $file, "rb" ) ) {
      die "Error: Cannot open '$file' - zlib error.\n";
    }
    return $fh;
  }
  elsif ( $file =~ m/\.(ann|snp|txt|tsv|csv)\z/xm ) {
    return path($file)->filehandle();
  }
  else {
    say STDERR "Error: Cannot open '$file', unrecognized file extension.";
    exit(1);
  }
}

