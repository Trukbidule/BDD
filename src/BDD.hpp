#pragma once

#include "truth_table.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>

/* These are just some hacks to hash std::pair (for the unique table).
 * You don't need to understand this part. */
namespace std
{
template<class T>
inline void hash_combine( size_t& seed, T const& v )
{
  seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template<>
struct hash<pair<uint32_t, uint32_t>>
{
  using argument_type = pair<uint32_t, uint32_t>;
  using result_type = size_t;
  result_type operator() ( argument_type const& in ) const
  {
    result_type seed = 0;
    hash_combine( seed, in.first );
    hash_combine( seed, in.second );
    return seed;
  }
};
}

class BDD
{
public:
  using index_t = uint32_t;
  /* Declaring `index_t` as an alias for an unsigned integer.
   * This is just for easier understanding of the code.
   * This datatype will be used for node indices. */

  using var_t = uint32_t;
  /* Similarly, declare `var_t` also as an alias for an unsigned integer.
   * This datatype will be used for representing variables. */

private:
  struct Node
  {
    var_t v; /* corresponding variable */
    index_t T; /* index of THEN child */
    index_t E; /* index of ELSE child */
  };

public:
  explicit BDD( uint32_t num_vars )
    : unique_table( num_vars ), num_invoke_not( 0u ), num_invoke_and( 0u ), num_invoke_or( 0u ), 
      num_invoke_xor( 0u ), num_invoke_ite( 0u )
  {
    nodes.emplace_back( Node({num_vars, 0, 0}) ); /* constant 0 */
    nodes.emplace_back( Node({num_vars, 1, 1}) ); /* constant 1 */
    /* `nodes` is initialized with two `Node`s representing the terminal (constant) nodes.
     * Their `v` is `num_vars` and their indices are 0 and 1.
     * (Note that the real variables range from 0 to `num_vars - 1`.)
     * Both of their children point to themselves, just for convenient representation.
     *
     * `unique_table` is initialized with `num_vars` empty maps. */
  }

  /**********************************************************/
  /***************** Basic Building Blocks ******************/
  /**********************************************************/

  uint32_t num_vars() const
  {
    return unique_table.size();
  }

  /* Get the (index of) constant node. */
  index_t constant( bool value ) const
  {
    return value ? 1 : 0;
  }

  /* Look up (if exist) or build (if not) the node with variable `var`,
   * THEN child `T`, and ELSE child `E`. */
  index_t unique( var_t var, index_t T, index_t E ){
    assert( var < num_vars() && "Variables range from 0 to `num_vars - 1`." );
    assert( T < nodes.size() && "Make sure the children exist." );
    assert( E < nodes.size() && "Make sure the children exist." );
    assert( nodes[T].v > var && "With static variable order, children can only be below the node." );
    assert( nodes[E].v > var && "With static variable order, children can only be below the node." );

    /* Reduction rule: Identical children */
    if ( T == E )
    {
      return T;
    }

    /* Look up in the unique table. */
    const auto it = unique_table[var].find( {T, E} );
    if ( it != unique_table[var].end() )
    {
      /* The required node already exists. Return it. */
      return it->second;
    }
    else
    {
      /* Create a new node and insert it to the unique table. */
      index_t const new_index = nodes.size();
      nodes.emplace_back( Node({var, T, E}) );
      unique_table[var][{T, E}] = new_index;
      return new_index;
    }
  }

  /* Return a node (represented with its index) of function F = x_var or F = ~x_var. */
  index_t literal( var_t var, bool complement = false )
  {
    return unique( var, constant( !complement ), constant( complement ) );
  }
  
  /********************* Ref Operations *********************/
  
  index_t ref( index_t f ){
      return f;
  }
  
  void deref(index_t f){
      
  }

  /**********************************************************/
  /********************* BDD Operations *********************/
  /**********************************************************/

  /* Compute ~f */
  index_t NOT( index_t f )
  {
    assert( f < nodes.size() && "Make sure f exists." );
    ++num_invoke_not;

    /* trivial cases */
    if ( f == constant( false ) )
    {
      return constant( true );
    }
    if ( f == constant( true ) )
    {
      return constant( false );
    }

    Node const& F = nodes[f];
    var_t x = F.v;
    index_t f0 = F.E, f1 = F.T;

    index_t const r0 = NOT( f0 );
    index_t const r1 = NOT( f1 );
    return unique( x, r1, r0 );
  }

  /* Compute f ^ g */
  index_t XOR( index_t f, index_t g )
  {
    assert( f < nodes.size() && "Make sure f exists." );
    assert( g < nodes.size() && "Make sure g exists." );
    ++num_invoke_xor;

    /* trivial cases */
    if ( f == g )
    {
      return constant( false );
    }
    if ( f == constant( false ) )
    {
      return g;
    }
    if ( g == constant( false ) )
    {
      return f;
    }
    if ( f == constant( true ) )
    {
      return NOT( g );
    }
    if ( g == constant( true ) )
    {
      return NOT( f );
    }
    if ( f == NOT( g ) )
    {
      return constant( true );
    }

    Node const& F = nodes[f];
    Node const& G = nodes[g];
    var_t x;
    index_t f0, f1, g0, g1;
    if ( F.v < G.v ) /* F is on top of G */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = g1 = g;
    }
    else if ( G.v < F.v ) /* G is on top of F */
    {
      x = G.v;
      f0 = f1 = f;
      g0 = G.E;
      g1 = G.T;
    }
    else /* F and G are at the same level */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = G.E;
      g1 = G.T;
    }

    index_t const r0 = XOR( f0, g0 );
    index_t const r1 = XOR( f1, g1 );
    return unique( x, r1, r0 );
  }

  /* Compute f & g */
  index_t AND( index_t f, index_t g )
  {
    assert( f < nodes.size() && "Make sure f exists." );
    assert( g < nodes.size() && "Make sure g exists." );
    ++num_invoke_and;

    /* trivial cases */
    if ( f == constant( false ) || g == constant( false ) )
    {
      return constant( false );
    }
    if ( f == constant( true ) )
    {
      return g;
    }
    if ( g == constant( true ) )
    {
      return f;
    }
    if ( f == g )
    {
      return f;
    }

    Node const& F = nodes[f];
    Node const& G = nodes[g];
    var_t x;
    index_t f0, f1, g0, g1;
    if ( F.v < G.v ) /* F is on top of G */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = g1 = g;
    }
    else if ( G.v < F.v ) /* G is on top of F */
    {
      x = G.v;
      f0 = f1 = f;
      g0 = G.E;
      g1 = G.T;
    }
    else /* F and G are at the same level */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = G.E;
      g1 = G.T;
    }

    index_t const r0 = AND( f0, g0 );
    index_t const r1 = AND( f1, g1 );
    return unique( x, r1, r0 );
  }

  /* Compute f | g */
  index_t OR( index_t f, index_t g )
  {
    assert( f < nodes.size() && "Make sure f exists." );
    assert( g < nodes.size() && "Make sure g exists." );
    ++num_invoke_or;

    /* trivial cases */
    if ( f == constant( true ) || g == constant( true ) )
    {
      return constant( true );
    }
    if ( f == constant( false ) )
    {
      return g;
    }
    if ( g == constant( false ) )
    {
      return f;
    }
    if ( f == g )
    {
      return f;
    }

    Node const& F = nodes[f];
    Node const& G = nodes[g];
    var_t x;
    index_t f0, f1, g0, g1;
    if ( F.v < G.v ) /* F is on top of G */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = g1 = g;
    }
    else if ( G.v < F.v ) /* G is on top of F */
    {
      x = G.v;
      f0 = f1 = f;
      g0 = G.E;
      g1 = G.T;
    }
    else /* F and G are at the same level */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      g0 = G.E;
      g1 = G.T;
    }

    index_t const r0 = OR( f0, g0 );
    index_t const r1 = OR( f1, g1 );
    return unique( x, r1, r0 );
  }

  /* Compute ITE(f, g, h), i.e., f ? g : h */
  index_t ITE( index_t f, index_t g, index_t h )
  {
    assert( f < nodes.size() && "Make sure f exists." );
    assert( g < nodes.size() && "Make sure g exists." );
    assert( h < nodes.size() && "Make sure h exists." );
    ++num_invoke_ite;

    /* trivial cases */
    if ( f == constant( true ) )
    {
      return g;
    }
    if ( f == constant( false ) )
    {
      return h;
    }
    if ( g == h )
    {
      return g;
    }

    Node const& F = nodes[f];
    Node const& G = nodes[g];
    Node const& H = nodes[h];
    var_t x;
    index_t f0, f1, g0, g1, h0, h1;
    if ( F.v <= G.v && F.v <= H.v ) /* F is not lower than both G and H */
    {
      x = F.v;
      f0 = F.E;
      f1 = F.T;
      if ( G.v == F.v )
      {
        g0 = G.E;
        g1 = G.T;
      }
      else
      {
        g0 = g1 = g;
      }
      if ( H.v == F.v )
      {
        h0 = H.E;
        h1 = H.T;
      }
      else
      {
        h0 = h1 = h;
      }
    }
    else /* F.v > min(G.v, H.v) */
    {
      f0 = f1 = f;
      if ( G.v < H.v )
      {
        x = G.v;
        g0 = G.E;
        g1 = G.T;
        h0 = h1 = h;
      }
      else if ( H.v < G.v )
      {
        x = H.v;
        g0 = g1 = g;
        h0 = H.E;
        h1 = H.T;
      }
      else /* G.v == H.v */
      {
        x = G.v;
        g0 = G.E;
        g1 = G.T;
        h0 = H.E;
        h1 = H.T;
      }
    }

    index_t const r0 = ITE( f0, g0, h0 );
    index_t const r1 = ITE( f1, g1, h1 );
    return unique( x, r1, r0 );
  }

  /**********************************************************/
  /***************** Printing and Evaluating ****************/
  /**********************************************************/

  /* Print the BDD rooted at node `f`. */
  void print( index_t f, std::ostream& os = std::cout ) const
  {
    for ( auto i = 0u; i < nodes[f].v; ++i )
    {
      os << "  ";
    }
    if ( f <= 1 )
    {
      os << "node " << f << ": constant " << f << std::endl;
    }
    else
    {
      os << "node " << f << ": var = " << nodes[f].v << ", T = " << nodes[f].T 
         << ", E = " << nodes[f].E << std::endl;
      for ( auto i = 0u; i < nodes[f].v; ++i )
      {
        os << "  ";
      }
      os << "> THEN branch" << std::endl;
      print( nodes[f].T, os );
      for ( auto i = 0u; i < nodes[f].v; ++i )
      {
        os << "  ";
      }
      os << "> ELSE branch" << std::endl;
      print( nodes[f].E, os );
    }
  }

  /* Get the truth table of the BDD rooted at node f. */
  Truth_Table get_tt( index_t f ) const
  {
    assert( f < nodes.size() && "Make sure f exists." );
    //assert( num_vars() <= 6 && "Truth_Table only supports functions of no greater than 6 variables." );
    if ( f == constant( false ) )
    {
      return Truth_Table( num_vars() );
    }
    else if ( f == constant( true ) )
    {
      return ~Truth_Table( num_vars() );
    }
    
    /* Shannon expansion: f = x f_x + x' f_x' */
    var_t const x = nodes[f].v;
    index_t const fx = nodes[f].T;
    index_t const fnx = nodes[f].E;
    //pb next line
    Truth_Table const tt_x = create_tt_nth_var( num_vars(), x );
    Truth_Table const tt_nx = create_tt_nth_var( num_vars(), x, false );
    return ( tt_x & get_tt( fx ) ) | ( tt_nx & get_tt( fnx ) );
  }

  /* Whether `f` is dead (having a reference count of 0). */
  bool is_dead( index_t f ) const
  {
    /* TODO */
    return false;
  }

  /* Get the number of living nodes in the whole package, excluding constants. */
  uint64_t num_nodes() const
  {
    uint64_t n = 0u;
    for ( auto i = 2u; i < nodes.size(); ++i )
    {
      if ( !is_dead( i ) )
      {
        ++n;
      }
    }
    return n;
  }

  /* Get the number of nodes in the sub-graph rooted at node f, excluding constants. */
  uint64_t num_nodes( index_t f ) const
  {
    assert( f < nodes.size() && "Make sure f exists." );

    if ( f == constant( false ) || f == constant( true ) )
    {
      return 0u;
    }

    std::vector<bool> visited( nodes.size(), false );
    visited[0] = true;
    visited[1] = true;

    return num_nodes_rec( f, visited );
  }

  uint64_t num_invoke() const
  {
    return num_invoke_not + num_invoke_and + num_invoke_or + num_invoke_xor + num_invoke_ite;
  }

private:
  /**********************************************************/
  /******************** Helper Functions ********************/
  /**********************************************************/

  uint64_t num_nodes_rec( index_t f, std::vector<bool>& visited ) const
  {
    assert( f < nodes.size() && "Make sure f exists." );
    

    uint64_t n = 0u;
    Node const& F = nodes[f];
    assert( F.T < nodes.size() && "Make sure the children exist." );
    assert( F.E < nodes.size() && "Make sure the children exist." );
    if ( !visited[F.T] )
    {
      n += num_nodes_rec( F.T, visited );
      visited[F.T] = true;
    }
    if ( !visited[F.E] )
    {
      n += num_nodes_rec( F.E, visited );
      visited[F.E] = true;
    }
    return n + 1u;
  }

private:
  std::vector<Node> nodes;
  std::vector<std::unordered_map<std::pair<index_t, index_t>, index_t>> unique_table;
  /* `unique_table` is a vector of `num_vars` maps storing the built nodes of each variable.
   * Each map maps from a pair of node indices (T, E) to a node index, if it exists.
   * See the implementation of `unique` for example usage. */

   /* Computed tables for each operation type. */
   std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_AND;
   std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_OR;
   std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_XOR;
   std::unordered_map<std::tuple<signal_t, signal_t, signal_t>, signal_t> computed_table_ITE;

  /* statistics */
  uint64_t num_invoke_not, num_invoke_and, num_invoke_or, num_invoke_xor, num_invoke_ite;
};
