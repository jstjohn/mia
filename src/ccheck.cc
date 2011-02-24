#include <map>
#include <string>
#include <utility>

#include <getopt.h>
#include <math.h>
#include <stdio.h>

extern "C" {
#include "map_align.h"
#include "mia.h"
#include "myers_align.h"
}

/*
 * Contamination Checker.  Outline:
 *
 * - read the human reference; this will contain ambiguity codes
 * - read maln file, including assembly and assembled reads
 * - align reference-consensus and assembly globally
 *   This uses Myers' O(nd) aligner, for it grasps ambiguity codes and
 *   runs fast enough for long, but similar sequences.
 * - find "diagnostic positions", positions where ass and ref differ
 * - for every "end" fragment: store it  and later join with its other
 *   half
 * - for every "full" fragment: if it crosses at least one diagnostic
 *   position, cut out that range from ref and align to it globally
 *   using the mia aligner
 * - for every position where the bases agree, classify it, then
 *   classify the fragment (conflicting, uninformative, contaminant,
 *   endogenous)
 * - produce a summary
 *
 * Notable features:
 * - uses substitution matrix from maln file
 * - operates sensibly on aDNA
 * - has sensible commandline and doesn't make too much noise in operation
 * - optionally considers only certain diagnostic positions
 *   (tranversions only and/or some region only)
 *
 * To be done:
 * - read multiple maln files if desired and sum up statistics
 *   [not needed right now and postponed]
 * - filter out some diagnostic positions depending on shape of sequence id
 *   [not needed right now and postponed]
 * - filter out some differences as not being diagnostic (esp. near
 *   homopolymers; see existing contamination checker)
 *   [low priority for now and postponed]
 *
 * Nice to have:
 * - check if small/capital letters and all ambiguity codes work
 *   everywhere sensibly, then make use of them
 * - somehow account for sequencing error, especially on Solexa, when
 *   calculating the likely contamination level
 * - remove linear scans where binary search methods would work
 */

void print_aln( const char* aln1, const char* aln2 )
{
	int p ;
	const char* a ;
	while( *aln1 && *aln2 ) {
		for( p = 0, a = aln1 ; *a && p != 72 ; ++p ) putchar( *a++ ) ;
		putchar( '\n' ) ;

		for( p = 0, a = aln2 ; *a && p != 72 ; ++p ) putchar( *a++ ) ;
		putchar( '\n' ) ;

		for( p = 0 ; *aln1 && *aln2 && p != 72 ; ++p )
			putchar( ( *aln1++ == *aln2++ ? '*' : ' ' ) ) ;
		putchar( '\n' ) ;
		putchar( '\n' ) ;
	}
}

// List of diagnostic positions: Coordinates are relative to assembly
// (we want to quickly know whether a fragment overlaps a DP).  We'll
// store the reference bases along with it.
typedef std::map< int, std::pair< char, char > > dp_list ;

// Everything that differs counts as diagnostic, unless it's a gap.  In
// principle, Ns could be diagnostic, too, even though in only one
// direction.  In practice, however, it turned out that Ns produce noise
// and little in the way of useable results.  So Ns don't count as
// diagnostic for now.
bool is_diagnostic( const char* aln1, const char* aln2 )
{
	return *aln1 != *aln2
				&& *aln1 != 'N' && *aln2 != 'N' 
				&& *aln1 != '-' && *aln2 != '-' ;
}

bool is_transversion( char a, char b )
{
	char u = a & ~32 ;
	char v = b & ~32 ;
	switch( u )
	{
		case 'A': return v != 'G' ;
		case 'C': return v != 'T' ;
		case 'G': return v != 'A' ;
		case 'T':
		case 'U': return v != 'C' ;
		default: return false ;
	}
}


dp_list mk_dp_list( const char* aln1, const char* aln2, bool transversions, int span_from, int span_to )
{
	dp_list l ;
	while( span_from != span_to && *aln1 && *aln2 )
	{
		if( is_diagnostic( aln1, aln2 ) && ( !transversions || is_transversion( *aln1, *aln2 )))
			l[span_from] = std::make_pair( *aln1, *aln2 ) ;
		if( *aln2 != '-' ) ++span_from ;
		++aln1 ;
		++aln2 ;
	}
	return l ;
}

std::pair< dp_list::const_iterator, dp_list::const_iterator >
overlapped_diagnostic_positions( const dp_list& l, const AlnSeqP s )
{
	dp_list::const_iterator left  = l.lower_bound( s->start ) ;
	dp_list::const_iterator right = l.lower_bound( s->end + 1 ) ;
	return std::make_pair( left, right ) ;
}

// XXX: linear scan --> O(n)
// This could be faster (O(log n)) if precompiled into some sort of index.
std::string lift_over( const char* aln1, const char* aln2, int s, int e ) 
{
	std::string r ;
	int p ;
	for( p = 0 ; p < e && *aln1 && *aln2 ; ++aln1, ++aln2 )
	{
		if( *aln1 != '-' && p >= s ) r.push_back( *aln1 ) ;
		if( *aln2 != '-' ) ++p ;
	}
	return r ;
}
 
bool consistent( bool adna, char x, char y )
{
	char x_ = x == 'G' ? 'R' : x == 'C' ? 'Y' : x ;
	return x == '-' || y == '-' || (char_to_bitmap( adna ? x_ : x ) & char_to_bitmap(y)) != 0 ;
}

enum whatsit { unknown, clean, dirt, conflict, nonsense, maxwhatsits } ;

const char *label[] = { "unclassified", "clean       ", "polluting   ", "conflicting ", "nonsensical " } ;

whatsit merge_whatsit( whatsit a, whatsit b )
{
	if( a == b ) return a ;
	if( a == unknown ) return b ;
	if( b == unknown ) return a ;
	if( a == nonsense || b == nonsense ) return nonsense ;
	return conflict ;
}

struct option longopts[] = {
	{ "reference", required_argument, 0, 'r' },
	{ "ancient", no_argument, 0, 'a' },
	{ "verbose", no_argument, 0, 'v' },
	{ "help",    no_argument, 0, 'h' },
	{ "transversions", no_argument, 0, 't' },
	{ "span", required_argument, 0, 's' },
	{ "maxd", required_argument, 0, 'd' },
	{ 0,0,0,0 }
} ;

void usage( const char* pname )
{
	fputs( "Usage: ", stdout ) ;
	fputs( pname, stdout ) ;
	fputs( " [-r <ref.fa>] [-a] [-v] <aln.maln> \n\n"
		"Reads a maln file and tries to quantify contained contamination.\n"
		"Options:\n"
		"  -r, --reference FILE     FASTA file with the likely contaminant\n"
		"  -a, --ancient            treat DNA as ancient (i.e. likely deaminated)\n"
		"  -t, --transversions      only transversions are diagnostic\n"
		"  -s, --span M-N           only look at range from M to N\n"
		"  -d, --maxd D             allow up to D differences between the references\n"
		"  -v, --verbose            increases verbosity level (can be repeated)\n"
		"  -h, --help               print this help message\n\n", stdout ) ;
}

int main( int argc, char * const argv[] )
{
	int summary[ maxwhatsits ] = {0} ;
	struct refseq hum_ref ;
	bool adna = false ;
	bool transversions = false ;
	int verbose = 0 ;
	int maxd = 1000 ;
	const char* ref_file = "mt311.fna" ;
	int span_from = 0, span_to = INT_MAX ;

	if( argc == 0 ) { usage( argv[0] ) ; return 0 ; }

	int opt ;
	do {
		opt = getopt_long( argc, argv, "r:avhts:d:", longopts, 0 ) ;
		switch( opt ) 
		{
			case 'r': 
				ref_file = optarg ;
				break ;
			case 'a':
				adna = true ;
				break ;
			case 'v':
				++verbose ;
				break ;
			case ':':
				puts( "missing option argument" ) ;
				break ;
			case '?':
				puts( "unknown option" ) ;
				break ;
			case 'h':
				usage( argv[0] ) ;
				return 0 ;
			case 't':
				transversions = true ;
				break ;
			case 's':
				sscanf( optarg, "%u-%u", &span_from, &span_to ) ;
				span_from-- ;
				break ;
			case 'd':
				maxd = atoi( optarg ) ;
				break ;
		}
	} while( opt != -1 ) ;

	if( optind == argc ) { usage( argv[0] ) ; return 0 ; }

	read_fasta_ref( &hum_ref, ref_file ) ;
	MapAlignmentP maln = read_ma( argv[optind] ) ;
	PSSMP submat = maln->fpsm ;

	char aln_con[ strlen(hum_ref.seq) + maxd + 1 ] ;
	char aln_ass[ strlen(maln->ref->seq) + maxd + 1 ] ;
	unsigned d = myers_diff( hum_ref.seq, myers_align_globally, maln->ref->seq, maxd, aln_con, aln_ass ) ;

	if( d == UINT_MAX ) { puts( "Couldn't align references (try to increase maxd)." ) ; return 1 ; }
	if( verbose >= 1 ) printf( "%d total differences between reference and assembly.\n", d ) ;
	if( verbose >= 6 ) print_aln( aln_con, aln_ass ) ;
	
	dp_list l = mk_dp_list( aln_con, aln_ass, transversions, span_from, span_to ) ;

	if( verbose >=1 ) 
	{
		int t = 0 ;
		for( dp_list::const_iterator i = l.begin() ; i != l.end() ; ++i )
			if( is_transversion( i->second.first, i->second.second ) ) ++t ;
		printf( "%d diagnostic positions, %d of which are transversions.\n", l.size(), t ) ;
	}
	if( verbose >= 3 ) 
	{
		dp_list::const_iterator i = l.begin() ;
		if( i != l.end() ) { printf( "<%d:%c,%c>", i->first, i->second.first, i->second.second ) ; ++i ; }
		for( ; i != l.end() ; ++i ) printf( ", <%d:%c,%c>", i->first, i->second.first, i->second.second ) ;
		putchar( '\n' ) ;
	}

	typedef std::map< std::string, std::pair< whatsit, int > > Bfrags ;
	Bfrags bfrags ;
	const AlnSeqP *s ;

	for( s = maln->AlnSeqArray ; s != maln->AlnSeqArray + maln->num_aln_seqs ; ++s )
	{
		whatsit klass = unknown ;
		int votes = 0 ;

		std::pair< dp_list::const_iterator, dp_list::const_iterator > p =
			overlapped_diagnostic_positions( l, *s ) ;
		if( p.first == p.second ) {
			if( verbose >= 3 ) {
				fputs( (*s)->id, stdout ) ;
				putchar( '/' ) ;
				putchar( (*s)->segment ) ;
				puts( ": no diagnostic positions" ) ;
			}
		}
		else
		{
			if( verbose >= 3 )
			{
				printf( "%s/%c: %d diagnostic positions", (*s)->id, (*s)->segment, std::distance( p.first, p.second ) ) ;
				if( verbose >= 4 ) 
				{
					putchar( ':' ) ; putchar( ' ' ) ;
					dp_list::const_iterator i = p.first ;
					if( i != p.second ) { printf( "<%d:%c,%c>", i->first, i->second.first, i->second.second ) ; ++i ; }
					for( ; i != p.second ; ++i ) printf( ", <%d:%c,%c>", i->first, i->second.first, i->second.second ) ;
				}
				printf( "\nrange:  %d..%d\n", (*s)->start, (*s)->end ) ;
			}

			std::string the_read ;
			for( char *nt = (*s)->seq, **ins = (*s)->ins ; *nt ; ++nt, ++ins )
			{
				if( *nt != '-' ) the_read.push_back( *nt ) ;
				if( *ins ) the_read.append( *ins ) ;
			}
			std::string the_ass( maln->ref->seq + (*s)->start, (*s)->end - (*s)->start + 1 ) ;
			std::string lifted = lift_over( aln_con, aln_ass, (*s)->start, (*s)->end + 1 ) ;

			if( verbose >= 5 )
			{
				printf( "raw read: %s\nlifted:   %s\nassembly: %s\n\n"
						"aln.read: %s\naln.assm: %s\nmatches:  %s",
						the_read.c_str(), lifted.c_str(), the_ass.c_str(), 
						(*s)->seq, the_ass.c_str() ) ;
				std::string::const_iterator b = the_ass.begin(), e = the_ass.end() ;
				const char* pc = (*s)->seq ;
				while( b != e && *pc ) putchar( *b++ == *pc++ ? '*' : ' ' ) ;
				putchar( '\n' ) ;
			}

			int size = std::max( lifted.size(), the_read.size() ) ;

			AlignmentP frag_aln = init_alignment( size, size, 0, 0 ) ;

			frag_aln->seq1 = lifted.c_str() ;
			frag_aln->seq2 = the_read.c_str() ;
			frag_aln->len1 = size ;
			frag_aln->len2 = size ;
			frag_aln->sg5 = 1 ;
			frag_aln->sg3 = 1 ;
			frag_aln->submat = submat ;
			pop_s1c_in_a( frag_aln ) ;
			pop_s2c_in_a( frag_aln ) ;
			dyn_prog( frag_aln ) ;

			pw_aln_frag pwaln ;
			max_sg_score( frag_aln ) ;			// ARGH!  This has a vital side-effect!!!
			find_align_begin( frag_aln ) ;  	//        And so has this...
			populate_pwaln_to_begin( frag_aln, &pwaln ) ;
			pwaln.start = frag_aln->abc;

			if( verbose >= 5 )
			{
				printf( "\naln.read: %s\naln.ref:  %s\nmatches:  ", pwaln.frag_seq, pwaln.ref_seq ) ;
				const char *pc = pwaln.frag_seq, *pd = pwaln.ref_seq ;
				while( *pc && *pd ) putchar( *pd++ == *pc++ ? '*' : ' ' ) ;
				putchar( '\n' ) ;
				putchar( '\n' ) ;
			}

			free_alignment( frag_aln ) ;

			char *paln1 = aln_con, *paln2 = aln_ass ;
			int ass_pos = 0 ;
			while( ass_pos != (*s)->start && *paln1 && *paln2 ) 
			{
				if( *paln2 != '-' ) ass_pos++ ;
				++paln1 ;
				++paln2 ;
			}

			std::string in_ref = lifted.substr( 0, pwaln.start ) ;
			in_ref.append( pwaln.ref_seq ) ;

			char *in_frag_v_ref = pwaln.frag_seq ;
			char *in_ass = maln->ref->seq + (*s)->start ;
			char *in_frag_v_ass = (*s)->seq ;

			if( *paln1 != in_ref[0] || *paln1 == '-' ) printf( "huh? (R+%d) %.10s %.10s\n", pwaln.start, paln1, in_ref.c_str() ) ;
			if( *paln2 != in_ass[0] && *paln2 != '-' ) printf( "huh? (A+%d) %.10s %.10s\n", pwaln.start, paln2, in_ass ) ;

			while( ass_pos != (*s)->end +1 && *paln1 && *paln2 && !in_ref.empty() && *in_ass && *in_frag_v_ass && *in_frag_v_ref )
			{
				if( is_diagnostic( paln1, paln2 ) ) {
					if( verbose >= 4 )
						printf( "diagnostic pos.: %d %c/%c %c/%c",
								ass_pos, in_ref[0], *in_frag_v_ref, *in_ass, *in_frag_v_ass ) ;
					if( *in_frag_v_ref != *in_frag_v_ass ) 
					{
						if( verbose >= 4 ) puts( "in disagreement." ) ;
					}
					else
					{
						bool maybe_clean = consistent( adna, *in_ass, *in_frag_v_ass ) ;
						bool maybe_dirt =  consistent( adna, in_ref[0], *in_frag_v_ref ) ;

						if( verbose >= 4 )
						{
							fputs( maybe_dirt  ? "" : "in", stdout ) ;
							fputs( " consistent/", stdout ) ;
							fputs( maybe_clean ? "" : "in", stdout ) ;
							fputs( " consistent\n", stdout ) ; 
						}

						if( maybe_clean && !maybe_dirt && klass == unknown ) klass = clean ;
						if( maybe_clean && !maybe_dirt && klass == dirt    ) klass = conflict ;
						if( !maybe_clean && maybe_dirt && klass == unknown ) klass = dirt ;
						if( !maybe_clean && maybe_dirt && klass == clean   ) klass = conflict ;
						if( !maybe_clean && !maybe_dirt )                    klass = nonsense ;
						if( maybe_dirt != maybe_clean ) votes++ ;
					}
				}

				if( *paln1 != '-' ) {
					do {
						in_ref=in_ref.substr(1) ;
						in_frag_v_ref++ ;
					} while( in_ref[0] == '-' ) ;
				}
				if( *paln2 != '-' ) {
					ass_pos++ ;
					do {
						in_ass++ ;
						in_frag_v_ass++ ;
					} while( *in_ass == '-' ) ;
				}
				++paln1 ;
				++paln2 ;
			}
		}

		Bfrags::const_iterator i = bfrags.find( (*s)->id ) ;

		switch( (*s)->segment )
		{
			case 'b':
				bfrags[ (*s)->id ] = std::make_pair( klass, votes ) ;
				break ;

			case 'f':
				if( i == bfrags.end() ) 
				{
					fputs( (*s)->id, stdout ) ;
					puts( "/f is missing its back." ) ;
				}
				else
				{
					votes += i->second.second ;
					klass = merge_whatsit( klass, i->second.first ) ;
				}
				
			case 'a':
				if( verbose >= 2 ) 
					printf( "%s is %s (%d votes)\n", (*s)->id, label[klass], votes ) ;
				if( verbose >= 3 ) putchar('\n') ;
				summary[klass]++ ;
				break ;

			default:
				fputs( "don't know how to handle fragment type ", stdout ) ;
				putchar( (*s)->segment ) ;
				putchar( '\n' ) ;
		}
	}

	puts( "\nSummary:" ) ;
	for( whatsit klass = unknown ; klass != maxwhatsits ; klass = (whatsit)( (int)klass +1 ) )
	{
		printf( "%s fragments: %d", label[klass], summary[klass] ) ;
		if( klass == dirt )
		{
			double z = 1.96 ; // this is Z_{0.975}, giving a 95% confidence interval (I hope...)
			double k = summary[dirt], n = k + summary[clean] ;
			double p_ = k / n ;
			double c = p_ + 0.5 * z * z / n ;
			double w = z * sqrt( p_ * (1-p_) / n + 0.25 * z * z / (n*n) ) ;
			double d = 1 + z * z / n ;

			printf( " (%.1f .. %.1f .. %.1f%)",
				100.0 * (c-w) / d,         		// lower bound of CI
				100.0 * p_,         			// ML estimate
				100.0 * (c+w) / d ) ;      		// upper bound of CI
		}
		putchar( '\n' ) ;
	}
	putchar( '\n' ) ;
}


