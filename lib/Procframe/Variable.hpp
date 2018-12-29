#pragma ident "$Id$"

/**
 * @file Variable.hpp
 * Class to define and handle 'descriptions' of GNSS variables.
 */

#ifndef VARIABLE_HPP
#define VARIABLE_HPP


#include "DataStructures.hpp"


namespace gpstk
{

    /** @addtogroup DataStructures */
    //@{


    /// Class to define and handle 'descriptions' of GNSS variables.
    class Variable
    {
    public:

        /// Default constructor for Variable
        Variable();


        /** Common constructor for Variable.
         *  By default, it is indexed by SourceID.
         *
         * @param type              TypeID of variable.
         * @param sourceIndexed     Whether this variable is SourceID-indexed
         *                          or not. By default, it IS SourceID-indexed.
         * @param satIndexed        Whether this variable is SatID-indexed
         *                          or not. By default, it is NOT SatID-indexed.
         * @param arcIndexed        Whether this variable is Arc-indexed
         *                          or not. By default, it is NOT Arc-indexed.
         * @param epochIndexed      Whether this variable is Epoch-indexed
         *                          or not. By default, it is NOT Epoch-indexed.
         * @param variance          Initial variance assigned to this variable.
         * @param currIndex         Current index in VariableSet.
         * @param prevIndex         Previous index in VariableSet.
         */
        Variable( const TypeID& type,
                  bool sourceIndexed    = true,
                  bool satIndexed       = false,
                  bool arcIndexed       = false,
                  bool epochIndexed     = false,
                  double variance         = 1e+10,
                  int  currIndex          = -1,
                  int  prevIndex          = -1 );


        /// Get variable type
        TypeID getType() const
        { return varType; };


        /** Set variable type
         *
         * @param type        New TypeID of variable.
         */
        Variable& setType(const TypeID& type)
        { varType = type; return (*this); };


        /// Get if this variable is SourceID-indexed
        bool getSourceIndexed() const
        { return isSourceIndexed; };


        /** Set if this variable is SourceID-indexed.
         *
         * @param sourceIndexed     Whether this variable is SourceID-indexed
         *                          or not. By default, it IS SourceID-indexed.
         */
        Variable& setSourceIndexed(bool sourceIndexed)
        { isSourceIndexed = sourceIndexed; return (*this); };


        /// Get if this variable is SatID-indexed.
        bool getSatIndexed() const
        { return isSatIndexed; };


        /** Set if this variable is SatID-indexed.
         *
         * @param satIndexed        Whether this variable is SatID-indexed
         *                          or not. By default, it is NOT SatID-indexed.
         */
        Variable& setSatIndexed(bool satIndexed)
        { isSatIndexed = satIndexed; return (*this); };


        /// Get if this variable is Type-indexed.
        bool getTypeIndexed() const
        { return isTypeIndexed; };


        /** Set if this variable is Type-indexed.
         *
         * @param typeIndexed       Whether this variable is Type-indexed
         *                          or not. By default, it IS Type-indexed.
         */
        Variable& setTypeIndexed(bool typeIndexed)
        { isTypeIndexed = typeIndexed; return (*this); };


        /// Get if this variable is Arc-indexed.
        bool getArcIndexed() const
        { return isArcIndexed; };


        /** Set if this variable is Arc-indexed.
         *
         * @param arcIndexed        Whether this variable is Arc-indexed
         *                          or not. By default, it IS Arc-indexed.
         */
        Variable& setArcIndexed(bool arcIndexed)
        { isArcIndexed = arcIndexed; return (*this); };


        /// Get if this variable is Epoch-indexed.
        bool getEpochIndexed() const
        { return isEpochIndexed; };


        /** Set if this variable is Epoch-indexed.
         *
         * @param epochIndexed      Whether this variable is Epoch-indexed
         *                          or not. By default, it IS Epoch-indexed.
         */
        Variable& setEpochIndexed(bool epochIndexed)
        { isEpochIndexed = epochIndexed; return (*this); };


        /// Get value of initial variance assigned to this variable.
        double getInitialVariance() const
        { return initialVariance; };


        /** Set value of initial variance assigned to this variable.
         *
         * @param variance      Initial variance assigned to this variable.
         */
        Variable& setInitialVariance(double variance)
        { initialVariance = variance; return (*this); };


        /// Get internal source this variable is assigned to (if any).
        SourceID getSource() const
        { return varSource; };


        /** Set internal source this variable is assigned to.
         *
         * @param source     Internal, specific SourceID of variable.
         */
        Variable& setSource(const SourceID& source)
        { varSource = source; return (*this); };


        /// Get internal satellite this variable is assigned to (if any).
        SatID getSatellite() const
        { return varSat; };


        /** Set internal satellite this variable is assigned to.
         *
         * @param satellite  Internal, specific SatID of variable.
         */
        Variable& setSatellite(const SatID& satellite)
        { varSat = satellite; return (*this); };


        /// Get internal arc this variable is assigned to (if any).
        int getArc() const
        { return varArc; };


        /** Set internal arc this variable is assigned to.
         *
         * @param arc  Internal, specific Arc of variable.
         */
        Variable& setArc(const int& arc)
        { varArc = arc; return (*this); };


        /// Get internal epoch this variable is assigned to (if any).
        CommonTime getEpoch() const
        { return varEpoch; };


        /** Set internal epoch this variable is assigned to.
         *
         * @param epoch  Internal, specific Epoch of variable.
         */
        Variable& setEpoch(const CommonTime& epoch)
        { varEpoch = epoch; return (*this); };


        /// Get index of current epoch
        int getCurrentIndex() const
        { return currentIndex; };


        /// Set index of current epoch
        void  setCurrentIndex(int index)
        { currentIndex = index; };


        /// Get index of previous epoch
        int getPreviousIndex() const
        { return previousIndex; };


        /// Set index of previous epoch
        void setPreviousIndex( int index)
        { previousIndex = index; };


        /// Equality operator
        virtual bool operator==(const Variable& right) const;


        /// This ordering is somewhat arbitrary, but is required to be able
        /// to use a Variable as an index to a std::map, or as part of a
        /// std::set.
        virtual bool operator<(const Variable& right) const;


        /// Inequality operator
        bool operator!=(const Variable& right) const
        { return !(operator==(right)); }


        /// Assignment operator
        virtual Variable& operator=(const Variable& right);


        /// SourceID object representing all sources : type(Unknown),
        /// sourceName("").
        static SourceID allSources;


        /// SourceID object representing "some" sources : type(Mixed),
        /// sourceName("").
        static SourceID someSources;


        /// SatID object representing no satellites:
        /// system(systemUnknown), id(-1).
        static SatID noSats;


        /// SatID object representing all satellites:
        /// system(systemMixed), id(-1).
        static SatID allSats;


        /// SatID object representing all satellites of GPS System:
        /// system(systemGPS), id(-1).
        static SatID allGPSSats;


        /// SatID object representing all satellites of Galileo System:
        /// system(systemGalileo), id(-1).
        static SatID allGalileoSats;


        /// SatID object representing all satellites of GLONASS System:
        /// system(systemGLONASS), id(-1).
        static SatID allGlonassSats;


        /// SatID object representing all satellites of BDS System:
        /// system(systemBDS), id(-1).
        static SatID allBDSSats;


        /// Destructor
        virtual ~Variable() {};


    private:

        /// Type of the variable
        TypeID varType;

        /// Source of the variable
        SourceID varSource;

        /// Sat of the variable
        SatID varSat;

        /// Arc of the variable
        int varArc;

        /// Epoch of the variable
        CommonTime varEpoch;


        bool isSourceIndexed;

        bool isSatIndexed;

        bool isTypeIndexed;

        bool isArcIndexed;

        bool isEpochIndexed;


        /// Value of initial variance assigned to this variable.
        double initialVariance;


        /// the index of this Variable in previous VariableSet.
        int previousIndex;


        /// the index of this Variable in current VariableSet.
        int currentIndex;


        /** Initializing function
         *
         * @param type        TypeID of variable.
         * @param pModel      Pointer to StochasticModel associated with
         *                    this variable. By default, it is a white
         *                    noise model.
         * @param variance    Initial variance assigned to this variable.
         */
        void Init( const TypeID& type,
                   double variance = 4e+14,
                   bool typeIndex  = true,
                   int currIndex   = -1,
                   int prevIndex   = -1);


    }; // End of class 'Variable'


    /// A structure used to store the coefficent information for a Variable.
    /// created by shjzhang
    struct Coefficient
    {
        /// Default constructor
        Coefficient( bool forceCoef   = false,
                     double coef      = 1.0  )
            : forceDefault(forceCoef), defaultCoefficient(coef)
        {};


        /// Boolean indicating if default coefficient is always used.
        bool forceDefault;


        /// Value of default coefficient assigned to this variable.
        double defaultCoefficient;


        /// Equality operator
        bool operator==(const Coefficient& right) const
        {
            return ( ( defaultCoefficient == right.defaultCoefficient ) &&
                     ( forceDefault == right.forceDefault ) );
        }


        /// Destructor
        virtual ~Coefficient() {};

    };


    /// Handy type definition

    typedef std::vector<Variable> VariableVector;

    typedef std::set<Variable> VariableSet;

    typedef std::list<Variable> VariableList;

    typedef std::map<Variable, double> VariableDataMap;

    typedef std::map<Variable, Coefficient> VarCoeffMap;


    namespace StringUtils
    {
        inline std::string asString(const Variable& v)
        {
            std::ostringstream oss;
            oss << v.getType() << "   "
                << v.getSource() << "   "
                << v.getSatellite() << "   "
                << v.getArc() << "   "
                << v.getEpoch() << "   "
                << v.getTypeIndexed() << " "
                << v.getSourceIndexed() << " "
                << v.getSatIndexed() << " "
                << v.getArcIndexed() << " "
                << v.getEpochIndexed();

            return oss.str();
        }
    }

    //@}

}  // End of namespace gpstk

#endif   // VARIABLE_HPP
