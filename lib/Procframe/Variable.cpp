#pragma ident "$Id$"

/**
 * @file Variable.cpp
 * Class to define and handle 'descriptions' of GNSS variables.
 */


#include "Variable.hpp"


namespace gpstk
{

    // SourceID object representing all sources : type(Unknown),
    // sourceName("").
    SourceID Variable::allSources;

    // SourceID object representing "some" sources : type(Mixed),
    // sourceName("").
    SourceID Variable::someSources(SourceID::Mixed, "");


    // SatID object representing no satellites:
    // system(systemUnknown), id(-1).
    SatID Variable::noSats( -1, SatID::systemUnknown );

    // SatID object representing all satellites:
    // system(systemMixed), id(-1).
    SatID Variable::allSats( -1, SatID::systemMixed );


    // SatID object representing all satellites of GPS System:
    // system(systemGPS), id(-1).
    SatID Variable::allGPSSats( -1, SatID::systemGPS );

    // SatID object representing all satellites of GLONASS System:
    // system(systemGLONASS), id(-1).
    SatID Variable::allGlonassSats( -1, SatID::systemGLONASS );

    // SatID object representing all satellites of Galileo System:
    // system(systemGalileo), id(-1).
    SatID Variable::allGalileoSats( -1, SatID::systemGalileo );

    // SatID object representing all satellites of BDS System:
    // system(systemBDS), id(-1).
    SatID Variable::allBDSSats(-1, SatID::systemBDS);


    // Default constructor for Variable
    Variable::Variable()
    {
        // Unknown variable type
        TypeID type;

        // Call Init method
        Init( type );

    }  // End of 'Variable::Variable()'



    /** Common constructor for Variable.
     *  By default, it is indexed by SourceID.
     *
     * @param type             TypeID of variable.
     * @param sourceIndexed    Whether this variable is SourceID-indexed
     *                         or not. By default, it IS SourceID-indexed.
     * @param satIndexed       Whether this variable is SatID-indexed
     *                         or not. By default, it is NOT.
     * @param variance         Initial variance assigned to this variable.
     * @param currIndex        Current index in VariableSet.
     * @param prevIndex        Previous index in VariableSet.
     */
    Variable::Variable( const TypeID& type,
                        bool sourceIndexed,
                        bool satIndexed,
                        bool arcIndexed,
                        bool epochIndexed,
                        double variance,
                        int  currIndex,
                        int  prevIndex )
    {
        bool typeIndex(true);

        // Call Init method
        Init(type, variance, typeIndex, currIndex, prevIndex);

        // This couple lines override settings by Init.
        isSourceIndexed = sourceIndexed;
        isSatIndexed = satIndexed;
        isArcIndexed = arcIndexed;
        isEpochIndexed = epochIndexed;

    }  // End of 'Variable::Variable()'



    /** Initializing function
     *
     * @param type        TypeID of variable.
     * @param variance    Initial variance assigned to this variable.
     * @param currIndex   Current index in VariableSet.
     * @param prevIndex   Previous index in VariableSet.
     */
    void Variable::Init( const TypeID& type,
                         double variance,
                         bool typeIndex,
                         int currIndex,
                         int prevIndex )
    {
        varType = type;

        // By default, it is source-indexed
        isSourceIndexed = true;

        // By default, it is not sat-indexed
        isSatIndexed = false;

        // By default, it is type-indexed
        isTypeIndexed = typeIndex;

        // By default, it is not arc-indexed
        isArcIndexed = false;

        // By default, it is not epoch-indexed
        isEpochIndexed = false;

        // Set initial variance
        initialVariance = variance;

        currentIndex = currIndex;

        previousIndex = prevIndex;

    }  // End of method 'Variable::Init()'



    // Equality operator
    bool Variable::operator==(const Variable& right) const
    {
        return ( ( varType == right.getType() )                     &&
                 ( isSourceIndexed == right.getSourceIndexed() )    &&
                 ( isSatIndexed == right.getSatIndexed() )          &&
                 ( isArcIndexed == right.getArcIndexed() )          &&
                 ( isEpochIndexed == right.getEpochIndexed() )      &&
//                 ( initialVariance == right.getInitialVariance() )  &&
                 ( varSource == right.getSource() )                 &&
                 ( varSat == right.getSatellite() )                 &&
                 ( varArc == right.getArc() )                       &&
                 ( varEpoch == right.getEpoch() )                   &&
                 ( isTypeIndexed == right.getTypeIndexed() ) );

    }  // End of 'Variable::operator=='



    // This ordering is somewhat arbitrary, but is required to be able
    // to use a Variable as an index to a std::map, or as part of a
    // std::set.
    bool Variable::operator<(const Variable& right) const
    {
        // Compare each field in turn
        if( varType == right.getType() )
        {
            // source-indexed
            if( isSourceIndexed == right.getSourceIndexed() )
            {
                // sat-indexed
                if( isSatIndexed == right.getSatIndexed() )
                {
                    // arc-indexed
                    if( isArcIndexed == right.getArcIndexed() )
                    {
                        // epoch-indexed
                        if( isEpochIndexed == right.getEpochIndexed() )
                        {
                            // initial variance
                            if( initialVariance == right.getInitialVariance() )
                            {
                                // source
                                if ( varSource == right.getSource() )
                                {
                                    // sat
                                    if( varSat == right.getSatellite() )
                                    {
                                        // arc
                                        if( varArc == right.getArc() )
                                        {
                                            // epoch
                                            if( varEpoch == right.getEpoch() )
                                            {
                                                return ( isTypeIndexed < right.getTypeIndexed() );
                                            }
                                            else
                                            {
                                                return ( varEpoch < right.getEpoch() );
                                            }
                                        }
                                        else
                                        {
                                            return ( varArc < right.varArc );
                                        }
                                    }
                                    else
                                    {
                                        return ( varSat < right.getSatellite() );
                                    }
                                }
                                else
                                {
                                    return ( varSource < right.getSource() );
                                }
                            }
                            else
                            {
                                return ( initialVariance < right.getInitialVariance() );
                            }
                        }
                        else
                        {
                            return (isEpochIndexed < right.getEpochIndexed() );
                        }
                    }
                    else
                    {
                        return ( isArcIndexed < right.getArcIndexed() );
                    }
                }
                else
                {
                    return ( isSatIndexed < right.getSatIndexed() );
                }
            }
            else
            {
                return ( isSourceIndexed < right.getSourceIndexed() );
            }
        }
        else
        {
            return ( varType < right.getType() );
        }

    }  // End of 'Variable::operator<'



    // Assignment operator
    Variable& Variable::operator=(const Variable& right)
    {

        // First check if these Variables are the same
        if ( this == &right ) return (*this);

        // If Variables are different, then set values of all fields
        setType( right.getType() );

        setSourceIndexed( right.getSourceIndexed() );

        setSatIndexed( right.getSatIndexed() );

        setArcIndexed( right.getArcIndexed() );

        setEpochIndexed( right.getEpochIndexed() );

        setInitialVariance( right.getInitialVariance() );

        setSource( right.getSource() );

        setSatellite( right.getSatellite() );

        setArc( right.getArc() );

        setEpoch( right.getEpoch() );

        setTypeIndexed(right.getTypeIndexed());

        setCurrentIndex(right.getCurrentIndex());

        setPreviousIndex(right.getPreviousIndex());

        return *this;

    }  // End of 'Variable::operator='

}  // End of namespace gpstk
