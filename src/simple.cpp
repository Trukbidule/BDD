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
    bdd.print( f );
    cout << tt << endl;
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

  {
    cout << "test 03: x0.x1.x2.x3.x4.x5.x6.x7 +x1.x8)" << endl;
    BDD bdd( 9 );
    auto const x0 = bdd.literal( 0 );
    auto const x1 = bdd.literal( 1 );
    auto const x2 = bdd.literal( 2 );
    auto const x3 = bdd.literal( 3 );
    auto const x4 = bdd.literal( 4 );
    auto const x5 = bdd.literal( 5 );
    auto const x6 = bdd.literal( 6 );
    auto const x7 = bdd.literal( 7 );
    auto const x8 = bdd.literal( 8 );
    
    auto const and01 = bdd.AND(x0, x1);
    auto const and23 = bdd.AND(x2, x3);
    auto const and45 = bdd.AND(x4, x5);
    auto const and67 = bdd.AND(x6, x7);
    auto const and18 = bdd.AND(x1, x8);
    auto const andTop1 = bdd.AND(and01, and23);
    auto const andTop2 = bdd.AND(and45, and67);
    auto const andTop3 = bdd.AND(andTop1, andTop2);
    
    auto const f = bdd.OR(andTop3, and18);
    
    auto const tt = bdd.get_tt( f );
    //bdd.print(f);
    cout << "Obtained: "<<endl;
    cout << tt << endl;
    string correction = "11001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001100110011001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
    cout << "Expected: "<<endl;
    cout << Truth_Table(correction)<<endl;
    passed &= check( tt, correction);
    //passed &= check( bdd.num_nodes( f ), 3 );
  }


  return passed ? 0 : 1;
}
