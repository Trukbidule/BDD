#pragma once

#include "truth_table.hpp"

#include <iostream>
#include <vector>
#include <unordered_map>
#include <functional>
#include <utility>
#include <tuple>
#include <string>
#include <vector>

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

class BDD{
public:
      using index_t = uint32_t;
      /* Declaring `index_t` as an alias for an unsigned integer.
       * This is just for easier understanding of the code.
       * This datatype will be used for node indices. */

      using var_t = uint32_t;
      /* Similarly, declare `var_t` also as an alias for an unsigned integer.
       * This datatype will be used for representing variables. */

       using signal_t = uint32_t;
       /* A signal represents an edge pointing to a node. 
       * The first 31 bits store the index of the node, 
       * and the last bit records whether the edge is complemented or not.
       * See below `make_signal`, `get_index` and `is_complemented`. */

private:
    struct Node{
        var_t v; /* corresponding variable */
        //children=signals: index of child+last bit=complemented
        signal_t Ts; /* index of THEN child */
        signal_t Es; /* index of ELSE child */
    };

    //create a signal from node index and complemented
    inline signal_t make_signal( index_t index, bool complement = false ) const{
        return complement ? ( index << 1 ) + 1 : index << 1;
        //+1 to indicate complemented
        //LSB=1 => bubble=complemented
    }

    //extract the node index from a signal
    inline index_t get_index( signal_t signal ) const{
        //std::cout<<"get_index: id: "<<( signal >> 1 )<<" size: "<<nodes.size() <<std::endl;
        assert( ( signal >> 1 ) < nodes.size() );
        return (signal >> 1);
    }
    
    //change only the complement of a signal
    inline signal_t set_complemented(signal_t signal, bool complement) const{
        return complement ? ( get_index(signal) << 1 ) + 1 : get_index(signal) << 1;
    }
    
    //toggle the complement of a signal
    inline signal_t toggle_complemented(signal_t signal) const{
        return is_complemented(signal) ? (signal&~0x1) : (signal|0x1);
        //return is_complemented(signal) ? (get_index(signal) << 1) : (( get_index(signal) << 1 ) + 1);
    }

    //extract complemented from a signal
    inline bool is_complemented( signal_t signal ) const{
        return (signal & 0x1);
    }
    
    //return the node from the signal (shortcut)
    inline Node get_node( signal_t signal ) const{
        //std::cout<<"get_node signal: "<<(signal>>1)<<std::endl;
        return nodes[get_index( signal )];
    }
  

public:
    explicit BDD( uint32_t num_vars )
    : unique_table( num_vars ), num_invoke_not( 0u ), num_invoke_and( 0u ), num_invoke_or( 0u ), 
    num_invoke_xor( 0u ), num_invoke_ite( 0u ){
        //old version:
        //nodes.emplace_back( Node({num_vars, 0, 0}) ); /* constant 0 */
        //nodes.emplace_back( Node({num_vars, 1, 1}) ); /* constant 1 */

        /* `nodes` is initialized with two `Node`s representing the terminal (constant) nodes.
        * Their `v` is `num_vars` and their indices are 0 and 1.
        * (Note that the real variables range from 0 to `num_vars - 1`.)
        * Both of their children point to themselves, just for convenient representation.
        *
        * `unique_table` is initialized with `num_vars` empty maps. */

        //new version: the constant 1 is at index 0, and its T and E point to itself
        nodes.emplace_back( Node({num_vars, 0, 0}) ); /* constant 1 */
    }

    /**********************************************************/
    /***************** Basic Building Blocks ******************/
    /**********************************************************/

    //getter
    uint32_t num_vars() const {
        return unique_table.size();
    }

    /* Get the (index of) constant node. */
    index_t constant_ind() const {
        //return value ? 1 : 0;
        return 0;//now only 1 at index 0
    }

    /* Get the (signal of) constant node. */
    signal_t constant_sig(bool value) const {
        return value ? make_signal(constant_ind(), 0) : make_signal(constant_ind(), 1);
    }
    
    //find the signals incomming on the node n
    std::vector<signal_t> get_signals_on(index_t n){
        std::vector<signal_t> sig;
        //scan the outgoing signal of all nodes, skip the constant at i=0
        for(index_t i=1; i<this->nodes.size(); i++){
            //if a signal from node i points to n: add it to sig (safety: avoid loop)
            if( (get_index(this->nodes[i].Ts) == n) && (get_index(this->nodes[i].Ts) != i) )
                sig.emplace_back(this->nodes[i].Ts);
            if( (get_index(this->nodes[i].Es) == n) && (get_index(this->nodes[i].Es) != i) )
                sig.emplace_back(this->nodes[i].Es);
        }
        return sig;
    }

    /* Look up (if exist) or build (if not) the node with variable `var`,
    * THEN child `T`, and ELSE child `E`. */
    signal_t unique( var_t var, signal_t Ts, signal_t Es ) {
        //std::cout<<"unique: var: "<<var<<" T: "<< (Ts>>1)<<" E: "<< (Es>>1)<<std::endl;
        
        //adapt from signal to index:
        index_t t = get_index(Ts);
        index_t e = get_index(Es);
        
        //std::cout<< "var:" <<std::to_string(var) << "num_var:" <<std::to_string(num_vars())<<std::endl;
        assert( var < num_vars() && "Variables range from 0 to `num_vars - 1`. var:");
        assert( t < nodes.size() && "Make sure the children exist." );
        assert( e < nodes.size() && "Make sure the children exist." );
        assert( nodes[t].v > var && "With static variable order, children can only be below the node." );
        assert( nodes[e].v > var && "With static variable order, children can only be below the node." );

        /* Reduction rule: Identical children */
        if ( Ts == Es ) {
            return Ts;
        }
        /* Reduction rule: remove the complement on T edge PART 1*/
        if( is_complemented(Ts) ){//if complement on Ts:
            toggle_complemented(Ts); //remove the complement of Ts
            toggle_complemented(Es); //toggle the complement of Es
            //toggle the complement of edges incomming on the node in part 2
        }


        /* Look up in the unique table. */
        const auto it = unique_table[var].find( {Ts, Es} );
        if ( it != unique_table[var].end() ){
            /* The required node already exists. Return it. */
            
            /* Reduction rule: remove the complement on T edge PART 2*/
            //if the node exist: toggle the complement of edges incomming onto it
            std::vector<signal_t> sig = get_signals_on(it->second);
            for(int i=0; i<sig.size(); i++){
                toggle_complemented( sig.at(i) );
            }
            
            return make_signal(it->second, false);
        } else {
            /* Create a new node and insert it to the unique table. */
            index_t const new_index = nodes.size();
            nodes.emplace_back( Node({var, Ts, Es}) );
            unique_table[var][{Ts, Es}] = new_index;
            return make_signal(new_index, false);
        }
    }

    /* Return a node (represented with its index) of function F = x_var or F = ~x_var. */
    signal_t literal( var_t var, bool complement = false ) {
        return unique( var, constant_sig( !complement ), constant_sig( complement ) );
    }
  
    /********************* Ref Operations *********************/
  
    index_t ref( index_t f ){
        //TODO PHASE 2
        return f;
    }

    void deref(index_t f){
        //TODO PHASE 2
    }

    /**********************************************************/
    /********************* BDD Operations *********************/
    /**********************************************************/

    /* Compute ~f */
    signal_t NOT(signal_t fs) {
        //signal_t sig = toggle_complemented(fs);
        return toggle_complemented(fs);
        //return unique( get_node(fs).v, toggle_complemented(get_node(fs).Ts), toggle_complemented(get_node(sig).Es) );
        //return unique( get_node(sig).v, get_node(sig).Ts, get_node(sig).Es );
        
        //adapt from signal to index:
        // index_t f = get_index(fs);
        // bool fc = is_complemented(fs);
        // 
        // assert( f < nodes.size() && "Make sure f exists." );
        // ++num_invoke_not;
        // 
        // /* trivial cases */
        // if ( f == constant() ) {//if cst: return opposit
        //     return make_signal(constant(), !fc);
        // }
        // 
        // Node const& F = nodes[f];
        // var_t x = F.v;
        // signal_t f0s = F.E, f1s = F.T;
        // 
        // signal_t const r0s = NOT(f0s);
        // signal_t const r1s = NOT(f1s);
        // return unique( x, r1s, r0s );
    }

    /* Compute f ^ g */
    signal_t XOR(signal_t fs, signal_t gs) {
        //std::cout<<"XOR f: "<<(fs>>1)<<" g: "<<(gs>>1)<<std::endl;
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        ++num_invoke_xor;

        /* trivial cases */
        if ( fs == gs ) {//x^x=0
            return constant_sig(false);
        }
        if ( fs == constant_sig(false) ) {//0^x=x
            return gs;
        }
        if ( gs == constant_sig(false) ) {//x^0=x
            return fs;
        }
        if ( fs == constant_sig(true) ) {//1^x=!x
            return NOT(gs);
        }
        if ( gs == constant_sig(true) ) {//x^1=!x
            return NOT(fs);
        }
        if ( f == g ) {//case f==!g, x^!x=1
            //the case fs==gs has already been evaluated before, so
            //if f==g => only the complement is different => f==!g
            return constant_sig(true);
        }

        Node const& F = nodes[f];
        Node const& G = nodes[g];
        var_t x;
        signal_t f0s, f1s, g0s, g1s;
        if( F.v < G.v ) {/* F is on top of G */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = g1s = gs;
        } else if( G.v < F.v ) {/* G is on top of F */
            x = G.v;
            f0s = f1s = fs;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        } else {/* F and G are at the same level */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        }

        signal_t const r0s = XOR(f0s, g0s);
        signal_t const r1s = XOR(f1s, g1s);
        return unique(x, r1s, r0s);
    }

    /* Compute f & g */
    signal_t AND(signal_t fs, signal_t gs) {
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        ++num_invoke_and;

        /* trivial cases */
        if( fs == constant_sig(false) || gs == constant_sig(false) ){//0&x=0
            return constant_sig(false);
        }
        if( fs == constant_sig(true) ){//1&x=x
            return gs;
        }
        if( gs == constant_sig(true) ){//x&1=x
            return fs;
        }
        if( fs == gs ){//x&x=x
            return fs;
        }
        if( f == g ){//case f==!g, x&!x=0
            //the case fs==gs has already been evaluated before, so
            //if f==g => only the complement is different => f==!g
            return constant_sig(false);
        }

        Node const& F = nodes[f];
        Node const& G = nodes[g];
        var_t x;
        signal_t f0s, f1s, g0s, g1s;
        if( F.v < G.v ) {/* F is on top of G */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = g1s = gs;
        } else if( G.v < F.v ) {/* G is on top of F */
            x = G.v;
            f0s = f1s = fs;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        } else {/* F and G are at the same level */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        }

        signal_t const r0s = AND( f0s, g0s );
        signal_t const r1s = AND( f1s, g1s );
        return unique(x, r1s, r0s);
    }

    /* Compute f | g */
    signal_t OR(signal_t fs, signal_t gs) {
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        ++num_invoke_or;

        /* trivial cases */
        if( fs == constant_sig(true) || gs == constant_sig(true) ){//1+x=1
            return constant_sig(true);
        }
        if( fs == constant_sig(false) ){//0+x=x
            return gs;
        }
        if( gs == constant_sig(false) ){//x+0=x
            return fs;
        }
        if( fs == gs ){//x+x=x
            return fs;
        }

        Node const& F = nodes[f];
        Node const& G = nodes[g];
        var_t x;
        signal_t f0s, f1s, g0s, g1s;
        if( F.v < G.v ) {/* F is on top of G */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = g1s = gs;
        } else if( G.v < F.v ) {/* G is on top of F */
            x = G.v;
            f0s = f1s = fs;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        } else {/* F and G are at the same level */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            g0s = gc ? toggle_complemented(G.Es) : G.Es;
            g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
        }

        signal_t const r0s = OR(f0s, g0s);
        signal_t const r1s = OR(f1s, g1s);
        return unique(x, r1s, r0s);
    }

    /* Compute ITE(f, g, h), i.e., f ? g : h */
    index_t ITE( signal_t fs, signal_t gs, signal_t hs ) {
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        index_t h = get_index(hs);
        bool hc = is_complemented(hs);
        
        //longer version: shannon exp: f = f.g + !f.h
        // signal_t sig = OR( AND(fs, gs), AND(NOT(fs), hs) );
        // return unique( get_node(sig).v, get_node(sig).Ts, get_node(sig).Es );
        
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        assert( h < nodes.size() && "Make sure h exists." );
        ++num_invoke_ite;
        
        /* trivial cases */
        if ( fs == constant_sig(true) ){//ITE(1, g, h)=g
            return gs;
        }
        if ( fs == constant_sig(false) ){//ITE(0, g, h)=0
            return hs;
        }
        if ( gs == hs ){//ITE(f, g, g)=g
            return gs;
        }
        
        Node const& F = nodes[f];
        Node const& G = nodes[g];
        Node const& H = nodes[h];
        var_t x;
        signal_t f0s, f1s, g0s, g1s, h0s, h1s;
        if ( F.v <= G.v && F.v <= H.v ) {/* F is not lower than both G and H */
            x = F.v;
            f0s = fc ? toggle_complemented(F.Es) : F.Es;
            f1s = fc ? toggle_complemented(F.Ts) : F.Ts;
            if ( G.v == F.v ){
                g0s = gc ? toggle_complemented(G.Es) : G.Es;
                g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
            } else {
                g0s = g1s = gs;
            }
            if ( H.v == F.v ) {
                h0s = hc ? toggle_complemented(H.Es) : H.Es;
                h1s = hc ? toggle_complemented(H.Ts) : H.Ts;
            } else {
                h0s = h1s = hs;
            }
        } else {/* F.v > min(G.v, H.v) */
            f0s = f1s = fs;
            if ( G.v < H.v ) {
                x = G.v;
                g0s = gc ? toggle_complemented(G.Es) : G.Es;
                g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
                h0s = h1s = hs;
            } else if ( H.v < G.v ) {
                x = H.v;
                g0s = g1s = gs;
                h0s = hc ? toggle_complemented(H.Es) : H.Es;
                h1s = hc ? toggle_complemented(H.Ts) : H.Ts;
            } else { /* G.v == H.v */
                x = G.v;
                g0s = gc ? toggle_complemented(G.Es) : G.Es;
                g1s = gc ? toggle_complemented(G.Ts) : G.Ts;
                h0s = hc ? toggle_complemented(H.Es) : H.Es;
                h1s = hc ? toggle_complemented(H.Ts) : H.Ts;
            }
        }
        
        index_t const r0s = ITE( f0s, g0s, h0s );
        index_t const r1s = ITE( f1s, g1s, h1s );
        return unique( x, r1s, r0s );
    }

    /**********************************************************/
    /***************** Printing and Evaluating ****************/
    /**********************************************************/

    /* Print the BDD rooted at node `f`. */
    void print( signal_t fs, std::ostream& os = std::cout ) const {
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        
        for ( auto i = 0u; i < nodes[f].v; ++i ){
            os << "  ";
        }
        if ( ( fs == constant_sig(true) ) || ( fs == constant_sig(false) ) ){//display cst, only 0
            os << "node " << f << ": constant 1" <<(fc ? " (c)":" (nc)")<< std::endl;
        } else {
            os << "node " << f <<(fc ? " (c)":" (nc)")<< ": var = " << nodes[f].v << ", T = " << get_index(nodes[f].Ts) << ", E = " << get_index(nodes[f].Es) << std::endl;
        
            for ( auto i = 0u; i < nodes[f].v; ++i ){
                os << "  ";
            }
            
            os << "> THEN branch" << std::endl;
            print( nodes[f].Ts, os );
            for ( auto i = 0u; i < nodes[f].v; ++i ){
                os << "  ";
            }
            os << "> ELSE branch" << std::endl;
            print( nodes[f].Es, os );
        }
    }

    /* Get the truth table of the BDD rooted at node f. */
    Truth_Table get_tt( signal_t fs ) const {
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        
        assert( f < nodes.size() && "Make sure f exists." );
        //assert( num_vars() <= 6 && "Truth_Table only supports functions of no greater than 6 variables." );
        
        if( fs == constant_sig(false) ) {
            return Truth_Table( num_vars() );//gives a TT full of 0
        }
        else if( fs == constant_sig(true) ) {
            return ~Truth_Table( num_vars() );//gives a TT full of 1
        }

        /* Shannon expansion: f = x f_x + x' f_x' */
        var_t const x = nodes[f].v;
        //if fs complemented in the first place: inverse the children TT
        signal_t const fxs = fc ? toggle_complemented(nodes[f].Ts) : nodes[f].Ts;
        signal_t const fnxs = fc ? toggle_complemented(nodes[f].Es) : nodes[f].Es;
        Truth_Table const tt_x = create_tt_nth_var( num_vars(), x );
        Truth_Table const tt_nx = create_tt_nth_var( num_vars(), x, false );
        return ( tt_x & get_tt( fxs ) ) | ( tt_nx & get_tt( fnxs ) );
    }

    /* Whether `f` is dead (having a reference count of 0). */
    bool is_dead( index_t f ) const{
        /* TODO PHASE 2*/
        return false;
    }

    /* Get the number of living nodes in the whole package, excluding constants. */
    uint64_t num_nodes() const {
        uint64_t n = 0u;
        for ( auto i = 2u; i < nodes.size(); ++i ){
            if ( !is_dead( i ) ){
                ++n;
            }
        }
        return n;
    }

    /* Get the number of nodes in the sub-graph rooted at node f, excluding constants. */
    uint64_t num_nodes(signal_t fs) const {
        //adapt from signal to index:
        index_t f = get_index(fs);
        
        assert( f < nodes.size() && "Make sure f exists." );

        //still working since constant(x) gives always 0
        if ( fs == constant_sig(false) || fs == constant_sig(true) ){
            return 0u;
        }

        std::vector<bool> visited( nodes.size(), false );
        visited[0] = true;
        visited[1] = true;

        return num_nodes_rec( fs, visited );
    }

    uint64_t num_invoke() const {
        return num_invoke_not + num_invoke_and + num_invoke_or + num_invoke_xor + num_invoke_ite;
    }

private:
    /**********************************************************/
    /******************** Helper Functions ********************/
    /**********************************************************/

    uint64_t num_nodes_rec( signal_t fs, std::vector<bool>& visited ) const {
        //adapt from signal to index:
        index_t f = get_index(fs);
        index_t t = get_index(nodes[f].Ts);
        index_t e = get_index(nodes[f].Es);
        
        assert( f < nodes.size() && "Make sure f exists." );

        uint64_t n = 0u;
        Node const& F = nodes[f];
        assert( t < nodes.size() && "Make sure the children exist." );
        assert( e < nodes.size() && "Make sure the children exist." );
        
        if ( !visited[t] ){
            n += num_nodes_rec( F.Ts, visited );
            visited[t] = true;
        }
        if ( !visited[e] ){
            n += num_nodes_rec( F.Es, visited );
            visited[e] = true;
        }
        return n + 1u;
    }

    private:
        std::vector<Node> nodes;
        std::vector<std::unordered_map<std::pair<index_t, index_t>, index_t>> unique_table;
        /* `unique_table` is a vector of `num_vars` maps storing the built nodes of each variable.
        * Each map maps from a pair of node indices (T, E) to a node index, if it exists.
        * See the implementation of `unique` for example usage. */

        // /* Computed tables for each operation type. */
        // std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_AND;
        // std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_OR;
        // std::unordered_map<std::tuple<signal_t, signal_t>, signal_t> computed_table_XOR;
        // std::unordered_map<std::tuple<signal_t, signal_t, signal_t>, signal_t> computed_table_ITE;

        /* statistics */
        uint64_t num_invoke_not, num_invoke_and, num_invoke_or, num_invoke_xor, num_invoke_ite;
};
