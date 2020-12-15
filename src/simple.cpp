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

bool check( uint64_t dd_size, uint64_t expected )
{
  cout << "  checking BDD size";
  if ( dd_size <= expected ) /* Using complemented edges can reduce, but will never increase BDD size. */
  {
    cout << "...passed." << endl;
    return true;
  }
  else
  {
    cout << "...failed. (expect " << expected << ", but get " << dd_size << " nodes)" << endl;
    return false;
  }
}

int main()
{
  bool passed = true;
  {
    cout << "test 00: x0 XOR x1" << endl;
    BDD bdd( 2 );
    auto const x0 = bdd.literal( 0 );
    auto const x1 = bdd.literal( 1 );
    auto const f = bdd.XOR( x0, x1 );
    auto const tt = bdd.get_tt( f );
    //bdd.print( f );
    //cout << tt << endl;
    passed &= check( tt, "0110" );
    passed &= check( bdd.num_nodes( f ), 3 );
  }

  {
    cout << "test 01: x0 AND x1" << endl;
    BDD bdd( 2 );
    auto const x0 = bdd.literal( 0 );
    auto const x1 = bdd.literal( 1 );
    auto const f = bdd.AND( x0, x1 );
    auto const tt = bdd.get_tt( f );
    //bdd.print( f );
    //cout << tt << endl;
    passed &= check( tt, "1000" );
    passed &= check( bdd.num_nodes( f ), 2 );
  }

  {
    cout << "test 02: ITE(x0, x1, x2)" << endl;
    BDD bdd( 3 );
    auto const x0 = bdd.literal( 0 );
    auto const x1 = bdd.literal( 1 );
    auto const x2 = bdd.literal( 2 );
    auto const f = bdd.ITE( x0, x1, x2 );
    auto const tt = bdd.get_tt( f );
    //bdd.print( f );
    //cout << tt << endl;
    passed &= check( tt, "11011000" );
    passed &= check( bdd.num_nodes( f ), 3 );
  }

  return passed ? 0 : 1;
}
