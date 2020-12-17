#include "BDD.hpp"
#include "truth_table.hpp"

#include <iostream>
#include <string>
#include <algorithm>

using namespace std;

bool check( Truth_Table const& tt, string const& ans )
{
  cout << "  checking function correctness";
  if ( tt == Truth_Table( ans ) )
  {
    cout << "...passed." << endl;
    return true;
  }
  else
  {
    cout << "...failed. (expect " << ans << ", but get " << tt << ")" << endl;
    return false;
  }
}

bool check( Truth_Table const& tt1, Truth_Table const& tt2 )
{
  cout << "  checking function correctness";
  if ( tt1 == tt2 )
  {
    cout << "...passed." << endl;
    return true;
  }
  else
  {
    cout << "...failed. (expect " << tt2 << ", but get " << tt1 << ")" << endl;
    return false;
  }
}

bool checkLE( uint64_t actual, uint64_t expected )
{
  if ( actual <= expected )
  {
    cout << "...passed." << endl;
    return true;
  }
  else
  {
    cout << "...failed. (expect <= " << expected << ", but get " << actual << ")" << endl;
    return false;
  }
}

bool checkEQ( uint64_t actual, uint64_t expected )
{
  if ( actual == expected )
  {
    cout << "...passed." << endl;
    return true;
  }
  else
  {
    cout << "...failed. (expect " << expected << ", but get " << actual << ")" << endl;
    return false;
  }
}

int main()
{
  bool passed = true;
  {
    cout << "test 00: large truth table";
    Truth_Table tt( "00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000" );

    if ( tt.n_var() == 7 )
    {
      cout << "...passed." << endl;
    }
    else
    {
      cout << "...failed." << endl;
      passed = false;
    }
  }

  {
    cout << "test 01: computed table" << endl;
    BDD bdd( 2 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const g = bdd.ref( bdd.AND( x0, x1 ) );
    auto const h = bdd.ref( bdd.AND( x0, x1 ) );
    bdd.deref( x0 ); bdd.deref( x1 );

    auto const f = bdd.ref( bdd.XOR( g, h ) );
    bdd.deref( g ); bdd.deref( h );
    auto const tt = bdd.get_tt( f );

    passed &= check( tt, "0000" );
    cout << "  checking number of computation";
    passed &= checkLE( bdd.num_invoke(), 5 );
  }

  {
    cout << "test 02: complemented edges" << endl;
    BDD bdd( 2 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const f = bdd.ref( bdd.XOR( x0, x1 ) );
    bdd.deref( x0 ); bdd.deref( x1 );
    auto const tt = bdd.get_tt( f );

    passed &= check( tt, "0110" );
    cout << "  checking BDD size (reachable nodes)";
    passed &= checkEQ( bdd.num_nodes( f ), 2 );

    cout << "test 03: reference count" << endl;
    cout << "  checking BDD size (living nodes)";
    passed &= checkEQ( bdd.num_nodes(), 2 );
  }

  {
    cout << "test 04: ITE(x2, x1, x0) AND ITE(x0, x2 AND NOT x1, x1 XOR x2)" << endl;
    BDD bdd( 3 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const x2 = bdd.ref( bdd.literal( 2 ) );

    auto const f1 = bdd.ref( bdd.ITE( x2, x1, x0 ) );

    auto const g = bdd.ref( bdd.AND( x2, bdd.NOT( x1 ) ) );
    auto const h = bdd.ref( bdd.XOR( x1, x2 ) );
    auto const f2 = bdd.ref( bdd.ITE( x0, g, h ) );
    bdd.deref( g ); bdd.deref( h );
    bdd.deref( x0 ); bdd.deref( x1 ); bdd.deref( x2 );

    auto const f = bdd.ref( bdd.AND( f1, f2 ) );
    bdd.deref( f1 ); bdd.deref( f2 );

    auto const tt = bdd.get_tt( f );
    passed &= check( tt, "00000000" );
    cout << "  checking BDD size (reachable nodes)";
    passed &= checkEQ( bdd.num_nodes( f ), 0 );
    cout << "  checking BDD size (living nodes)";
    passed &= checkEQ( bdd.num_nodes(), 0 );
  }

  {
    cout << "test 05: ITE(x2 AND x3, x1 AND NOT x0, NOT x2 AND NOT x4)" << endl;
    BDD bdd( 5 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const x2 = bdd.ref( bdd.literal( 2 ) );
    auto const x3 = bdd.ref( bdd.literal( 3 ) );
    auto const x4 = bdd.ref( bdd.literal( 4 ) );
    auto const f1 = bdd.ref( bdd.AND( x2, x3 ) );
    auto const f2 = bdd.ref( bdd.AND( x1, bdd.NOT( x0 ) ) );
    auto const f3 = bdd.ref( bdd.AND( bdd.NOT( x2 ), bdd.NOT( x4 ) ) );
    bdd.deref( x0 ); bdd.deref( x1 ); bdd.deref( x2 ); bdd.deref( x3 ); bdd.deref( x4 );

    auto const f = bdd.ref( bdd.ITE( f1, f2, f3 ) );
    bdd.deref( f1 ); bdd.deref( f2 ); bdd.deref( f3 );

    auto const tt = bdd.get_tt( f );
    passed &= check( tt, "01000000000000000100111100001111" );
    cout << "  checking BDD size (reachable nodes)";
    passed &= checkEQ( bdd.num_nodes( f ), 6 );
    cout << "  checking BDD size (living nodes)";
    passed &= checkEQ( bdd.num_nodes(), 6 );
  }

  {
    cout << "test 06: more than 6 variables & multiple POs" << endl;
    BDD bdd( 10 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const x2 = bdd.ref( bdd.literal( 2 ) );
    auto const x3 = bdd.ref( bdd.literal( 3 ) );
    auto const x4 = bdd.ref( bdd.literal( 4 ) );
    auto const x5 = bdd.ref( bdd.literal( 5 ) );
    auto const x6 = bdd.ref( bdd.literal( 6 ) );
    auto const x7 = bdd.ref( bdd.literal( 7 ) );
    auto const x8 = bdd.ref( bdd.literal( 8 ) );
    auto const x9 = bdd.ref( bdd.literal( 9 ) );

    auto const f1 = bdd.ref( bdd.OR( x0, x9 ) );

    auto const g1 = bdd.ref( bdd.AND( x6, bdd.NOT( x4 ) ) );
    auto const g2 = bdd.ref( bdd.AND( x4, bdd.NOT( x6 ) ) );
    auto const f2 = bdd.ref( bdd.OR( g1, g2 ) );
    bdd.deref( g1 ); bdd.deref( g2 );

    auto const f3 = bdd.ref( bdd.ITE( x6, bdd.NOT( x2 ), bdd.NOT( x6 ) ) );

    bdd.deref( x0 ); bdd.deref( x1 ); bdd.deref( x2 ); bdd.deref( x3 ); bdd.deref( x4 );
    bdd.deref( x5 ); bdd.deref( x6 ); bdd.deref( x7 ); bdd.deref( x8 ); bdd.deref( x9 );

    auto const tt1 = bdd.get_tt( f1 );
    passed &= check( tt1, create_tt_nth_var( 10, 0 ) | create_tt_nth_var( 10, 9 ) );
    auto const tt2 = bdd.get_tt( f2 );
    passed &= check( tt2, create_tt_nth_var( 10, 4 ) ^ create_tt_nth_var( 10, 6 ) );
    auto const tt3 = bdd.get_tt( f3 );
    passed &= check( tt3, ~create_tt_nth_var( 10, 2 ) | ~create_tt_nth_var( 10, 6 ) );

    cout << "  checking BDD size (reachable nodes) of f1";
    passed &= checkEQ( bdd.num_nodes( f1 ), 2 );
    cout << "  checking BDD size (reachable nodes) of f2";
    passed &= checkEQ( bdd.num_nodes( f2 ), 2 );
    cout << "  checking BDD size (reachable nodes) of f3";
    passed &= checkEQ( bdd.num_nodes( f3 ), 2 );
    cout << "  checking BDD size (living nodes)";
    passed &= checkEQ( bdd.num_nodes(), 5 );
  }

  {
    cout << "test 07: computed table for XOR" << endl;
    BDD bdd( 4 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const x2 = bdd.ref( bdd.literal( 2 ) );
    auto const x3 = bdd.ref( bdd.literal( 3 ) );

    auto const g1 = bdd.ref( bdd.XOR( x2, x3 ) );
    auto const g2 = bdd.ref( bdd.XOR( x1, g1 ) );
    auto const g3 = bdd.ref( bdd.XOR( x0, g2 ) );
    bdd.deref( g1 ); bdd.deref( g2 );

    auto const h1 = bdd.ref( bdd.XOR( x3, x2 ) );
    auto const h2 = bdd.ref( bdd.XOR( x0, x1 ) );
    auto const h3 = bdd.ref( bdd.XOR( h1, h2 ) );
    bdd.deref( h1 ); bdd.deref( h2 );

    bdd.deref( x0 ); bdd.deref( x1 ); bdd.deref( x2 ); bdd.deref( x3 );
    
    auto const f = bdd.ref( bdd.XOR( g3, h3 ) );
    bdd.deref( g3 ); bdd.deref( h3 );

    auto const tt = bdd.get_tt( f );
    passed &= check( tt, "0000000000000000" );
    cout << "  checking number of computation";
    passed &= checkLE( bdd.num_invoke(), 20 );
  }

  {
    cout << "test 08: computed table for ITE" << endl;
    BDD bdd( 3 );
    auto const x0 = bdd.ref( bdd.literal( 0 ) );
    auto const x1 = bdd.ref( bdd.literal( 1 ) );
    auto const x2 = bdd.ref( bdd.literal( 2 ) );

    auto const f1 = bdd.ref( bdd.ITE( x1, x2, x0 ) );
    auto const f2 = bdd.ref( bdd.ITE( bdd.NOT( x1 ), x0, x2 ) );
    bdd.deref( x0 ); bdd.deref( x1 ); bdd.deref( x2 );

    auto const tt1 = bdd.get_tt( f1 );
    passed &= check( tt1, "11100010" );
    auto const tt2 = bdd.get_tt( f2 );
    passed &= check( tt2, "11100010" );
    cout << "  checking number of computation";
    passed &= checkLE( bdd.num_invoke(), 10 );
  }

  if ( passed )
  {
    cout << endl << "All tests passed, congrats!" << endl;
  }

  return passed ? 0 : 1;
}
