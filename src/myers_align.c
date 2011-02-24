#include "myers_align.h"

#include <limits.h>
#include <string.h>

int min( int a, int b ) { return a < b ? a : b ; }
int max( int a, int b ) { return a < b ? b : a ; }
int max3( int a, int b, int c ) { return a < b ? max( b, c ) : max( a, c ) ; }

inline int match( char a, char b ) { return (char_to_bitmap(a) & char_to_bitmap(b)) != 0 ; }

// [*blech*, this looks and feels like FORTRAN.]
unsigned myers_diff( const char *seq_a, enum myers_align_mode mode, const char* seq_b, unsigned maxd, char *bt_a, char *bt_b ) 
{
	int len_a = strlen( seq_a ), len_b = strlen( seq_b ) ;
	if( maxd > len_a + len_b ) maxd = len_a + len_b ;
	int v_[ 2 * maxd + 1 ][maxd] ; // Hmmm...
	int (*v)[maxd] = v_ + maxd ;
	
	int d, dd, k, x, y ;
	for( d = 0 ; d != maxd ; ++d )		// D-paths in order of increasing D
	{
		for( k = max(-d,-len_a) ; k <= min(d,len_b) ; ++k ) // diagonals
		{
			if( d == 0 )         x = 0 ;
			else if( k == -d   ) x = v[ k+1 ] [d-1] ;
			else if( k == -d+1 ) x = max( v[ k ][d-1]+1, v[ k+1 ][d-1] ) ;
			else if( k ==  d   ) x = v[ k-1 ][d-1]+1 ;
			else if( k ==  d-1 ) x = max( v[ k ][d-1]+1, v[ k-1 ][d-1]+1 ) ;
			else                 x = max3( v[ k-1 ][d-1]+1, v[ k ][d-1]+1, v[ k+1 ][d-1] ) ;

			y = x-k ;
			while( x < len_b && y < len_a && match( seq_b[x], seq_a[y] ) ) ++x, ++y ;
			v[k][d] = x ;

			if(
					(mode == myers_align_is_prefix || y == len_a) &&
					(mode == myers_align_has_prefix || x == len_b) )
			{
				char *out_a = bt_a + len_a + d +1 ;
				char *out_b = bt_b + len_b + d +1 ;
				*--out_a = 0 ;
				*--out_a = 0 ;
				for( dd = d ; dd != 0 ; )
				{
					if( k != -dd && k != dd && x == v[ k ][dd-1]+1 )
					{
						--dd ;
						--x ;
						--y ;
						*--out_b = seq_b[x] ;
						*--out_a = seq_a[y] ;
					}
					else if( k > -dd+1 && x == v[ k-1 ][dd-1]+1 )
					{
						--x ;
						--k ;
						--dd ;
						*--out_b = seq_b[x] ;
						*--out_a = '-' ;
					}
					else if( k < dd-1 && x == v[ k+1 ][dd-1] )
					{
						++k ;
						--y ;
						--dd ;
						*--out_b = '-' ;
						*--out_a = seq_a[y] ;
					}
					else // this better had been a match...
					{
						--x ;
						--y ;
						*--out_b = seq_b[x] ;
						*--out_a = seq_a[y] ;
					}
				}
				while( x > 0 )
				{
					--x ;
					*--out_b = seq_b[x] ;
					*--out_a = seq_a[x] ;
				}
				memmove( bt_a, out_a, bt_a + len_a + d + 1 - out_a ) ;
				memmove( bt_b, out_b, bt_b + len_b + d + 1 - out_b ) ;
				return d ;
			}
		}
	}
	return UINT_MAX ;
}

