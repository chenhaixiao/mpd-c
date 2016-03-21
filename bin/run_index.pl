#!/usr/bin/env perl
# Name:           run_index.pl
# Date Created:   Mon Mar 21 13:59:19 2016
# Date Modified:  Mon Mar 21 13:59:19 2016
# By:             TS Wingo
#
# Description:

use 5.10.0;
use warnings;
use strict;

use Getopt::Long;
use Path::Tiny;
use Data::Dump qw/ dump /;

# variables
my ( $dir_name, $out_ext );
my $ext        = "fa";
my $indexDepth = 14;

# get options
die "Usage: $0 [-e <fasta extension>] [-i index depth] -d <dir_name> -o <out_ext>\n"
  unless GetOptions(
  'o|out=s'   => \$out_ext,
  'e|ext=s'   => \$ext,
  'd|dir=s'   => \$dir_name,
  'i|index=n' => \$indexDepth,
  )
  and $dir_name
  and $indexDepth
  and $out_ext
  and $ext;

my ( $filesAref, $sizesAref ) = FastaList( $dir_name, $ext );
WriteIn( $out_ext, $filesAref, $sizesAref );
WriteSh($out_ext);

sub WriteIn {
  my $out_ext   = shift;
  my $filesAref = shift;
  my $sizesAref = shift;

  my @sSizes = sort { $b <=> $a } @$sizesAref;
  my $maxSize = shift @sSizes;

  my $fh = path("$out_ext.in")->filehandle(">");
  say {$fh} join "\n", ( "d", "in.index_genome.in", scalar @$filesAref );
  say {$fh} join "\n", @$filesAref;
  say {$fh} "$out_ext.d$indexDepth";
  say {$fh} ( $maxSize + 1 );
  say {$fh} $indexDepth;
}

sub WriteSh {
  my $out_ext = shift;
  my $fh      = path("$out_ext.sh")->filehandle(">");
  say {$fh} qq{#!/bin/sh
./index_genome < $out_ext.in};
}

sub FastaList {
  my $dir = shift;
  my $ext = shift;

  my ( @fastas, @sizes );

  my $pt    = path($dir);
  my @files = path($dir)->children(qr{$ext\z});
  my @chrs  = ( 1 .. 26, "M", "X", "Y", "Un" );
  my %files = map { $_->basename() => $_ } @files;

  for my $chr (@chrs) {
    my $f = sprintf( "chr%s.fa", $chr );
    if ( exists $files{$f} ) {
      push @fastas, $files{$f}->stringify;
      push @sizes, ( -s $files{$f} );
      delete $files{$f};
    }
  }

  for my $chr (@chrs) {
    for my $file ( sort keys %files ) {
      if ( $file =~ m/\Achr$chr/ ) {
        push @fastas, $files{$file}->stringify;
        push @sizes, ( -s $files{$file} );
        delete $files{$file};
      }
    }
  }
  return ( \@fastas, \@sizes );
}
