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

//Template definitions, for the hash function (of unique_table and computed_table)
namespace std{
    template<class T>
        inline void hash_combine( size_t& seed, T const& v ){
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }

    template<>//for the unique_table
        struct hash<pair<uint32_t, uint32_t>>{
            using argument_type = pair<uint32_t, uint32_t>;
            using result_type = size_t;
            result_type operator() ( argument_type const& in ) const{
                result_type seed = 0;
                hash_combine( seed, in.first );
                hash_combine( seed, in.second );
                return seed;
            }
        };
    
    template<>//for the computed_tables
        struct hash<tuple<uint32_t, uint32_t>>{
            using argument_type = tuple<uint32_t, uint32_t>;
            using result_type = size_t;
            result_type operator() ( argument_type const& in ) const{
                result_type seed = 0;
                hash_combine( seed, get<0>(in) );
                hash_combine( seed, get<1>(in) );
                return seed;
            }
        };
        
    template<>//for the ITE computed_table
        struct hash<tuple<uint32_t, uint32_t, uint32_t>>{
            using argument_type = tuple<uint32_t, uint32_t, uint32_t>;
            using result_type = size_t;
            result_type operator() ( argument_type const& in ) const{
                result_type seed = 0;
                hash_combine( seed, get<0>(in) );
                hash_combine( seed, get<1>(in) );
                hash_combine( seed, get<2>(in) );
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
        assert( ( signal >> 1 ) < nodes.size() );
        return (signal >> 1);
    }
    
    //toggle the complement of a signal
    inline signal_t toggle_complemented(signal_t signal) const{
        return is_complemented(signal) ? (signal&~0x1) : (signal|0x1);
    }
    
    //change only the complement of a given signal
    inline signal_t set_complemented(signal_t signal, bool comp) const{
        return comp ? (signal | 0x1) : (signal & ~(0x1));
    }

    //extract complemented from a signal
    inline bool is_complemented( signal_t signal ) const{
        return (signal & 0x1);
    }
    
    //return the node from the signal (shortcut)
    inline Node get_node( signal_t signal ) const{
        return nodes[get_index( signal )];
    }
  

public:
    explicit BDD( uint32_t num_vars )
    : unique_table( num_vars ), num_invoke_not( 0u ), num_invoke_and( 0u ), num_invoke_or( 0u ), 
    num_invoke_xor( 0u ), num_invoke_ite( 0u ){

        /* `nodes` is initialized with two `Node`s representing the terminal (constant) nodes.
        * Their `v` is `num_vars` and their indices are 0 and 1.
        * (Note that the real variables range from 0 to `num_vars - 1`.)
        * Both of their children point to themselves, just for convenient representation.
        *
        * `unique_table` is initialized with `num_vars` empty maps. */

        //new version: the constant 1 is at index 0, and its T and E point to itself
        nodes.emplace_back( Node({num_vars, 0, 0}) ); /* constant 1 */
        //create the constant reference count(required to ensure same index for the others)
        ref_count.emplace_back(0);
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
        return 0;//now the only constant is 1, placed at index 0, with the var=num_var
    }

    /* Get the (signal of) constant node. */
    signal_t constant_sig(bool value) const {
        return value ? make_signal(constant_ind(), 0) : make_signal(constant_ind(), 1);
    }
    
    //find the signals incomming on the node n
    // std::vector<signal_t*> get_signals_on(index_t n){
    //     std::vector<signal_t*> sig;
    //     //scan the outgoing signal of all nodes, skip the constant at i=0
    //     for(index_t i=1; i<this->nodes.size(); i++){
    //         //if a signal from node i points to n: add it to sig (safety: avoid loop)
    //         if( (get_index(this->nodes[i].Ts) == n) && (get_index(this->nodes[i].Ts) != i) )
    //             sig.emplace_back( &(this->nodes[i].Ts) );
    //         if( (get_index(this->nodes[i].Es) == n) && (get_index(this->nodes[i].Es) != i) )
    //             sig.emplace_back( &(this->nodes[i].Es) );
    //     }
    //     return sig;
    // }
    
    //find the node having as a child the node n
    // std::vector<index_t> get_nodes_to(index_t n){
    //     std::vector<index_t> ret;
    //     //scan the outgoing signal of all nodes, skip the constant at i=0
    //     for(index_t i=1; i<this->nodes.size(); i++){
    //         //if a node has as child n, add it to the list
    //         if(get_index(this->nodes[i].Ts) == n)
    //             ret.emplace_back(i);
    //         else if(get_index(this->nodes[i].Es) == n)
    //             ret.emplace_back(i);
    //     }
    //     return ret;
    // }

    /* Look up (if exist) or build (if not) the node with variable `var`,
    * THEN child `T`, and ELSE child `E`. */
    signal_t unique( var_t var, signal_t Ts, signal_t Es ) {
        //std::cout<<"\t unique: var: "<<var<<" T: "<< (Ts>>1)<<" (v: "<<get_node(Ts).v<<") "<<" E: "<< (Es>>1)<<" (v: "<<get_node(Es).v<<") "<<std::endl;
        
        //adapt from signal to index:
        index_t t = get_index(Ts);
        index_t e = get_index(Es);
        
        assert( var < num_vars() && "Variables range from 0 to `num_vars - 1`. var:");
        assert( t < nodes.size() && "Make sure the children exist." );
        assert( e < nodes.size() && "Make sure the children exist." );
        
        assert( nodes[t].v > var && "With static variable order, children can only be below the node." );
        assert( nodes[e].v > var && "With static variable order, children can only be below the node." );
        
        //final values
        signal_t newTs = Ts;
        signal_t newEs = Es;
        bool comp = false;//by default: not complemented
        
        /* Reduction rule: remove the complement on T edge */
        if( is_complemented(Ts) ){//if complement on Ts:
            newTs = toggle_complemented(Ts); //remove the complement of Ts
            newEs = toggle_complemented(Es); //toggle the complement of Es
            //toggle the complement of edges incomming on the node in part 2
            comp = true;
        }
        
        /* Reduction rule: Identical children*/
        if(Ts == Es ) {//if identical children
            return make_signal(get_index(Ts), comp); //skip the node
        }

        /* Look up in the unique table. */
        const auto it = unique_table[var].find( {newTs, newEs} );
        if ( it != unique_table[var].end()){
            /* If the required node already exists. Return it. */
            return make_signal(it->second, comp);
        } else {
            /* Create a new node and insert it to the unique table. */
            index_t const new_index = nodes.size();
            nodes.emplace_back( Node({var, newTs, newEs}) );
            ref_count.emplace_back(0);//creat a reference for this node
            
            unique_table[var][{newTs, newEs}] = new_index;
            return make_signal(new_index, comp);
        }
    }

    /* Return a node (represented with its index) of function F = x_var or F = ~x_var. */
    signal_t literal( var_t var, bool complement = false ) {
        return unique( var, constant_sig( !complement ), constant_sig( complement ) );
    }
  
    /********************* Ref Operations *********************/
  
    //mark a signal and all its children as "in use"
    signal_t ref( signal_t fs ){
        //assert( get_index(fs) < num_vars());
        if( get_index(fs) != constant_ind() ){//if not the final node
            ref_count.at(get_index(fs))++;//increase its ref_count
            //do the same to its children
            ref(get_node(fs).Ts);
            ref(get_node(fs).Es);
        }
        return fs;
    }

    //remove the "in use" mark of a signal and its children
    void deref(signal_t fs){
        if( get_index(fs) != constant_ind() ){//if not the final node
            //decrease its ref_count if not already zero
            if(ref_count.at(get_index(fs)) > 0)
                ref_count.at(get_index(fs))--;
            //do the same to its children
            deref(get_node(fs).Ts);
            deref(get_node(fs).Es);
        }
    }

    /**********************************************************/
    /********************* BDD Operations *********************/
    /**********************************************************/

    /* Compute ~f */
    signal_t NOT(signal_t fs) {
        return toggle_complemented(fs);
    }

    /* Compute f ^ g */
    signal_t XOR(signal_t fs, signal_t gs) {
        ++num_invoke_xor;
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        //safety check
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        
        //check if already computed
        auto it = computed_table_XOR.find( std::make_tuple(fs, gs) );
        if( it != computed_table_XOR.end() ){//if already computed, return it
            return it->second;
        }
        it = computed_table_XOR.find( std::make_tuple(gs, fs) );
        if( it != computed_table_XOR.end() ){//if already computed, return it
            return it->second;
        }

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
        if ( fs == NOT(gs) ) {//case f==!g, x^!x=1
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
        signal_t ret = unique(x, r1s, r0s);
        //add to computed table
        computed_table_XOR.emplace(std::make_pair(std::make_tuple(fs, gs), ret ) );
        return ret;
    }

    /* Compute f & g */
    signal_t AND(signal_t fs, signal_t gs) {
        ++num_invoke_and;
        
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        //safety check
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        
        //check if already computed
        auto it = computed_table_AND.find( std::make_tuple(fs, gs) );
        if( it != computed_table_AND.end() ){//if already computed, return it
            return it->second;
        }
        it = computed_table_AND.find( std::make_tuple(gs, fs) );
        if( it != computed_table_AND.end() ){//if already computed, return it
            return it->second;
        }

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
        signal_t ret = unique(x, r1s, r0s);
        //add to computed table
        computed_table_AND.emplace(std::make_pair(std::make_tuple(fs, gs), ret ) );
        return ret;
    }

    /* Compute f | g */
    signal_t OR(signal_t fs, signal_t gs) {
        ++num_invoke_or;
        
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        
        //safety check
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        
        //check if already computed
        auto it = computed_table_OR.find( std::make_tuple(fs, gs) );
        if( it != computed_table_OR.end() ){//if already computed, return it
            return it->second;
        }
        it = computed_table_OR.find( std::make_tuple(gs, fs) );
        if( it != computed_table_OR.end() ){//if already computed, return it
            return it->second;
        }
        
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
        signal_t ret = unique(x, r1s, r0s);
        //add to computed table
        computed_table_OR.emplace(std::make_pair(std::make_tuple(fs, gs), ret) );
        return ret;
    }

    /* Compute ITE(f, g, h), i.e., f ? g : h */
    signal_t ITE( signal_t fs, signal_t gs, signal_t hs ) {
        ++num_invoke_ite;
        
        //adapt from signal to index:
        index_t f = get_index(fs);
        bool fc = is_complemented(fs);
        index_t g = get_index(gs);
        bool gc = is_complemented(gs);
        index_t h = get_index(hs);
        bool hc = is_complemented(hs);
        
        assert( f < nodes.size() && "Make sure f exists." );
        assert( g < nodes.size() && "Make sure g exists." );
        assert( h < nodes.size() && "Make sure h exists." );
        
        //check if already computed
        auto it = computed_table_ITE.find( std::make_tuple(fs, gs, hs) );
        if( it != computed_table_ITE.end() ){//if already computed, return it
            return it->second;
        }
        it = computed_table_ITE.find( std::make_tuple(NOT(fs), hs, gs) );
        if( it != computed_table_ITE.end() ){//if already computed, return it
            return it->second;
        }
        
        //test version, slower: shannon exp: f = f.g + !f.h
        // signal_t sig = OR( AND(fs, gs), AND(NOT(fs), hs) );
        // return unique( get_node(sig).v, get_node(sig).Ts, get_node(sig).Es );

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
        
        signal_t const r0s = ITE( f0s, g0s, h0s );
        signal_t const r1s = ITE( f1s, g1s, h1s );
        signal_t ret = unique( x, r1s, r0s );
        //add to computed table
        computed_table_ITE.emplace(std::make_pair(std::make_tuple(fs, gs, hs), ret) );
        return ret;
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
        //selection of the variables
        Truth_Table const tt_x = create_tt_nth_var( num_vars(), x );
        Truth_Table const tt_nx = create_tt_nth_var( num_vars(), x, false );
        //build the final TT
        return ( tt_x & get_tt( fxs ) ) | ( tt_nx & get_tt( fnxs ) );
    }

    /* Whether `f` is dead (having a reference count of 0). */
    bool is_dead( index_t f ) const{
        if(ref_count.at(f) != 0)
            return false;
        else
            return true;
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
        std::vector<uint32_t> ref_count;
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
