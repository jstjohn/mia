#!/usr/bin/env perl

use Getopt::Std;
use vars qw( $opt_r $opt_L $opt_t $opt_H $opt_m $opt_k $opt_w $opt_M );
use strict;
my $whmia_exec = `which mia`;
$whmia_exec =~ s/\s+$//;
my $whma_exec  = `which ma`;
$whma_exec =~ s/\s+$//;

my $ITER_NUM = 1;
my $MAX_ITER = 30;
&init();

my ( $cmd, $status );

system( "touch CONS.0.fa" );
$cmd = "$opt_w -L $opt_L -t $opt_t -r $opt_r -H $opt_H -m $opt_m -k $opt_k";
$status = system( $cmd );
if ( $status ) {
    print STDERR ( "$status from\n: $cmd\n" );
}
$cmd = ( "$opt_M -M $opt_m.1 -f 5 > CONS.$ITER_NUM.fa" );
$status = system( $cmd );
if ( $status ) {
    print STDERR ( "$status from\n $cmd\n" );
}

while( !&convergence( $opt_m ) &&
       ($ITER_NUM < $MAX_ITER) ) {

    # Reassemble
    $cmd = "$opt_w -r CONS.$ITER_NUM.fa -t $opt_t -L $opt_L -H $opt_H -m $opt_m -k $opt_k";
    $status = system( $cmd );
    if ( $status ) {
	print STDERR ( "$status from\n: $cmd\n" );
    }

    # Call new consensus
    $ITER_NUM++;
    $cmd = ( "$opt_M -M $opt_m.1 -f 5 > CONS.$ITER_NUM.fa" );
    $status = system( $cmd );
    if ( $status ) {
	print STDERR ( "$status from\n $cmd\n" );
    }
}

sub convergence {
    my $fn_root = shift;
    my $diff_lines = 0;
    my $diff_cmd;
    my $last_iter_num = $ITER_NUM - 1;

    # Construct diff command
    $diff_cmd = "diff CONS.$ITER_NUM.fa CONS.$last_iter_num.fa";

    # Is this latest consensus the same as before?
    open( DL, "$diff_cmd|" ) or die( "Can't do $diff_cmd" );
    
    while( <DL> ) {
	$diff_lines++;
    }

    close( DL );

    if ( $diff_lines == 0 ) {
	return 1;
    }
    else {
	return 0;
    }
}

sub init {
    getopts( 'r:L:t:H:m:k:w:M:' );
    unless( -f $opt_r &&
	    -f $opt_L ) {
	print( "iterate_whmia_list_D.pl\n" );
	print( "Iteratively run mia or wormhole-mia with -D option\n" );
	print( "OPTIONS:\n" );
	print( "-r [beginning reference sequence]\n" );
	print( "-L [input sequences file list to assemble]\n" );
	print( "-t [threads]; default = 1\n" );
	print( "-H [-H option for mia]\n" );
	print( "-m [output maln filename]\n" );
	print( "-k [-k option for mia]\n" );
	print( "-w [wormhole mia executable; default = $whmia_exec\n" );
	print( "-M [wormhole ma executable; default = $whma_exec\n" );
	exit( 0 );
    }

    unless( defined( $opt_w ) ) {
	$opt_w = $whmia_exec;
    }

    unless( defined( $opt_M ) ) {
	$opt_M = $whma_exec;
    }
}
    
