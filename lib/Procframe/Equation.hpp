#pragma ident "$Id$"

/**
 * @file Equation.hpp
 * GNSS Data Structure to define and handle 'descriptions' of GNSS equations.
 */

#ifndef GPSTK_EQUATION_HPP
#define GPSTK_EQUATION_HPP

//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
//  any later version.
//
//  The GPSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GPSTk; if not, write to the Free Software Foundation,
//  Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009
//
//============================================================================


#include "Variable.hpp"
#include <iterator>


namespace gpstk
{

    /** @addtogroup DataStructures */
    //@{


    /// Defines a header containing basic equation data
    struct equationHeader
    {

        /// Source this equation is related to
        SourceID equationSource;


        /// Satellite this equation is related to
        SatID equationSat;


        /// Epoch this equation is related to
        CommonTime equationEpoch;


        /** In case this variable is associated to SOME specific
         *  sources ("Variable::someSources" in "equationSource"),
         *  then the corresponding SourceID set is stored here.
         */
        std::set<SourceID> equationSourceSet;


        /** In case this variable is associated to SOME specific
         *  satellites if the size of the variable > 0.
         */
        std::set<SatID> equationSatSet;


        /// Independent term
        Variable indTerm;


        /** Constant weight associated to this equation. This is a relative
         *  value that compares with the other Equations. It is 1.0 by
         *  default.
         */
        double constWeight;


        /** Map to store the values related to type. It is retrieved from
         *  the gnssDataMap, and used to compute the state transition matrix
         *  and the process noise matrix.
         */
        typeValueMap typeValueData;


        /// Default constructor
        equationHeader()
            : equationSource(Variable::allSources),
              equationSat(Variable::allSats),
              equationEpoch(CommonTime::BEGINNING_OF_TIME),
              constWeight(1.0)
        {};


        /** Explicit constructor
         *
         * @param source    Source this equation is related to.
         * @param sat       Satellite this equation is related to.
         * @param epoch     Epoch this equation is related to.
         * @param var       Variable representing the independent term.
         * @param weight    Constant weight associated to this equation.
         */
        equationHeader( const SourceID& source,
                        const SatID& sat,
                        const CommonTime& epoch,
                        const Variable var,
                        const double& weight )
         : equationSource(source),
           equationSat(sat),
           equationEpoch(epoch),
           indTerm(var),
           constWeight(weight)
        {};


        /** Explicit constructor from a Variable
         *
         * @param var       Variable representing the independent term.
         */
        equationHeader(const Variable& var)
            : equationSource(Variable::allSources),
              equationSat(Variable::allSats),
              equationEpoch(CommonTime::BEGINNING_OF_TIME),
              indTerm(var),
              constWeight(1.0)
        {};


        /// Assignment operator
        virtual equationHeader& operator=(const equationHeader& right);


        /** Assignment operator from a Variable
         *
         * @param var       Variable representing the independent term.
         */
        virtual equationHeader& operator=(const Variable& var)
        { indTerm = var; return (*this); };


        /// Destructor
        virtual ~equationHeader() {};

    }; // End of struct 'equationHeader'



    /** GNSS Data Structure to define and handle 'descriptions' of GNSS
     *  equations.
     */
    struct Equation : gnssData<equationHeader, VariableDataMap>
    {

        /// Default constructor.
        Equation();


        /** Common constructor. It defines an Equation from its header. You
         *  must later use other methods to input the variables.
         *
         * @param head     Data structure describing the Equation header.
         */
        Equation( const equationHeader& head )
        { header = head; };


        /** Common constructor. It defines an Equation from its independent
         *  term. You must later use other methods to input the variables.
         *
         * @param indep     Variable object describing the independent term.
         */
        Equation( const Variable& indep );


        /** Common constructor. It defines an Equation from its independent
         *  term. You must later use other methods to input the variables.
         *
         * @param var     TypeID object describing the independent term.
         */
        Equation( const TypeID& type );


        /** Common constructor. It takes a simple gnssEquationDefinition
         *  object and creates a  more complex Equation object.
         *
         * @param gnssEq  gnssEquationDefinition object.
         *
         * A gnssEquationDefinition object defines equations as a simple list
         * of TypeID's: The independent term (usually the prefit residual)
         * type in the header, and the variables' types in the body (or
         * 'unknowns').
         *
         * The resulting Equation object will honor this simple structure,
         * assigning white noise models to all variables, as well as declaring
         * them source-specific and satellite-independent (or 'UN-specific').
         *
         * The former is suitable for simple GNSS data processing strategies
         * like the SPS C1-based positioning, where the variables are
         * TypeID::dx, TypeID::dy, TypeID::dz and TypeID::cdt.
         *
         */
        Equation( const gnssEquationDefinition& gnssEq );


        /// Return the independent term of this equation.
        virtual Variable getIndependentTerm() const
        { return header.indTerm; };


        /** Set the independent term of this Equation.
         *
         * @param var     Variable object describing the independent term.
         */
        virtual Equation& setIndependentTerm(const Variable& var)
        { header = var; return (*this); };


        /// Get the value of the constant weight associated to this equation.
        virtual double getWeight() const
        { return header.constWeight; };


        /** Set the value of the constant weight associated to this equation.
         *
         * @param cweight    Value of constant weight.
         */
        virtual Equation& setWeight(const double& cweight)
        { header.constWeight = cweight; return (*this); };


        /** Add a variable (unknown) to this Equation.
         *
         * @param var     Variable object to be added to the unknowns.
         */
        virtual Equation& addVariable(const Variable& var,
                                      const double coef = 1.0)
        {
            // Insert new variable and coefData into the body
            body.insert( std::make_pair(var, coef) );

            return (*this);
        };


        /** Add a variable (unknown) to this Equation.
         *
         * @param var     Variable object to be added to the unknowns.
         */
//        virtual Equation& addVariable( const Variable& var,
//                                       const bool forceCoef = false,
//                                       const double coef    = 1.0 )
//        {
//            // Coefficient object
//            Coefficient coefData(forceCoef, coef);
//
//            // Insert new variable and coefData into the body
//            body.insert( std::make_pair(var, coefData) );
//
//            return (*this);
//        };


        /** Add a variable (unknown) to this Equation.
         *
         * @param type             TypeID of variable.
         * @param pModel           Pointer to StochasticModel associated with
         *                         this variable. By default, it is a white
         *                         noise model.
         * @param sourceIndexed    Whether this variable is SourceID-indexed
         *                         or not. By default, it IS SourceID-indexed.
         * @param satIndexed       Whether this variable is SatID-indexed
         *                         or not. By default, it is NOT SatID-indexed.
         * @param variance         Initial variance assigned to this variable.
         * @param coef             Default coefficient assigned.
         */
        virtual Equation& addVariable( const TypeID& type,
                                       bool sourceIndexed       = true,
                                       bool satIndexed          = false,
                                       bool arcIndexed          = false,
                                       bool epochIndexed        = false,
                                       double variance          = 4.0e14,
                                       double coef              = 1.0 );


        /** Remove a variable (unknown) from this Equation.
         *
         * @param var     Variable object to be romoved from the unknowns.
         */
        virtual Equation& removeVariable(const Variable& var)
        { body.erase(var); return (*this); };


        /** Remove ALL variables (unknowns) from this Equation.
         *
         * @warning This method does NOT clear the Equation's independent
         *          term. You MUST take care of it yourself (use method
         *          'setIndependentTerm', for instance).
         */
        virtual Equation& clear()
        { body.clear(); return (*this); };


        /// Get equation SourceID.
        SourceID getEquationSource() const
        { return header.equationSource; };

        /// Get equation SatID.
        SatID getEquationSat() const
        { return header.equationSat; };


        /** Get SourceID set. This is only meaningful if "equationSource" in
         *  header is set to "Variable::someSources".
         */
        std::set<SourceID> getSourceSet() const
        { return header.equationSourceSet; };


         /** Get SatID set. This is only meaningful if size of 'equationSatSet'
          * is not zero.
          */
        std::set<SatID> getSatSet() const
        { return header.equationSatSet; };


        CommonTime getEquationEpoch() const
        { return header.equationEpoch; };


        /** Add a source to SourceID set. This is only meaningful if
         * "equationSource" in header is set to "Variable::someSources".
         */
        Equation& addSource2Set( const SourceID& source )
        { header.equationSourceSet.insert(source); return (*this); };


        /** Add a sat to SatID set. This is only meaningful if size of
         * 'equationSatSet' in header is not zero.
         */
        Equation& addSat2Set( const SatID& sat )
        { header.equationSatSet.insert(sat); return (*this); };


        /** Clear SourceID set. This is only meaningful if "equationSource"
         *  in header is set to "Variable::someSources".
         */
        Equation& clearSourceSet()
        { header.equationSourceSet.clear(); return (*this); };


        /** Clear SatID set. This is only meaningful if "equationSource"
         *  in header is set to "Variable::someSources".
         */
        Equation& clearSatSet()
        { header.equationSatSet.clear(); return (*this); };


        /// This ordering is somewhat arbitrary, but is required to be able
        /// to use an Equation as an index to a std::map, or as part of a
        /// std::set.
        virtual bool operator<(const Equation& right) const
        { return (header.indTerm < right.header.indTerm); };


        /// Destructor
        virtual ~Equation() {};


    }; // End of struct 'Equation'


    /// Handy type definition

    typedef std::list<Equation> EquationList;

    typedef std::map<SourceID, EquationList> SourceEquationMap;


    namespace StringUtils
    {
        inline std::string asString(const VariableDataMap& vmap)
        {
            std::ostringstream oss;
            for( VariableDataMap::const_iterator it = vmap.begin();
                 it != vmap.end();
                 ++it )
            {
                oss << (*it).first.getType() << "   "
                    << (*it).first.getSource() << "   "
                    << (*it).first.getSatellite() << "   "
                    << (*it).first.getEpoch() << "   "
                    << (*it).first.getTypeIndexed() << " "
                    << (*it).first.getSourceIndexed() << " "
                    << (*it).first.getSatIndexed() << " "
                    << (*it).first.getArcIndexed() << " "
                    << (*it).first.getEpochIndexed()
                    << std::endl;
            }

            return oss.str();
        }
    }

    //@}

}  // End of namespace gpstk

#endif   // GPSTK_EQUATION_HPP
