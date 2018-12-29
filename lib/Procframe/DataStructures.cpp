#pragma ident "$Id$"

/**
 * @file DataStructures.cpp
 * Set of several data structures to be used by other GPSTk classes.
 */

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
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Dagoberto Salazar - gAGE ( http://www.gage.es ). 2007, 2008, 2009, 2011
//
//============================================================================
//
//  Revision
//
//  2012/06/19   Add 'getGnssRinex(SatID& sat)' method for class gnssDataMap.
//               The new method will return the first gnssRinex data which
//               contain the given satellite. ( 2012/06/19 )
//
//  2014/03/16   add the sourceNumber for 'source' in header.
//               f.header.source.sourceNumber= roh.markerNumber;
//
//  2014/07/15   Define a temporary 'source' to initialize the value of
//               the source in header; (shjzhang)
//
//============================================================================


#include "DataStructures.hpp"

using namespace gpstk::StringUtils;
using namespace std;

namespace gpstk
{


      ////// typeValueMap //////


      // Return a TypeIDSet with all the data types present in
      // this object.
    TypeIDSet typeValueMap::getTypeID() const
    {

        TypeIDSet typeSet;

        for( typeValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            typeSet.insert( (*pos).first );
        }

        return typeSet;

    }  // End of method 'typeValueMap::getTypeID()'



      // Return a typeValueMap with only this type of data.
      // @param type Type of value to be extracted.
    typeValueMap typeValueMap::extractTypeID(const TypeID& type) const
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return extractTypeID(typeSet);

    }  // End of method 'typeValueMap::extractTypeID()'



      // Return a typeValueMap with only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data to
      //                be extracted.
    typeValueMap typeValueMap::extractTypeID(const TypeIDSet& typeSet) const
    {

        typeValueMap tvMap;

        for( TypeIDSet::const_iterator pos = typeSet.begin();
             pos != typeSet.end();
             ++pos )
        {

            typeValueMap::const_iterator itObs( (*this).find(*pos) );
            if( itObs != (*this).end() )
            {
                tvMap[ (*itObs).first ] = (*itObs).second;
            }
        }

        return tvMap;

    }  // End of method 'typeValueMap::extractTypeID()'



      // Modifies this object, keeping only this type of data.
      // @param type Type of value to be kept.
    typeValueMap& typeValueMap::keepOnlyTypeID(const TypeID& type)
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return (keepOnlyTypeID(typeSet));

    }  // End of method 'typeValueMap::keepOnlyTypeID()'



      // Modifies this object, keeping only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    typeValueMap& typeValueMap::keepOnlyTypeID(const TypeIDSet& typeSet)
    {

        typeValueMap tvMap( (*this).extractTypeID(typeSet) );
        (*this) = tvMap;

        return (*this);

    }  // End of method 'typeValueMap::keepOnlyTypeID()'



      // Modifies this object, removing these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    typeValueMap& typeValueMap::removeTypeID(const TypeIDSet& typeSet)
    {

        for( TypeIDSet::const_iterator pos = typeSet.begin();
             pos != typeSet.end();
             ++pos )
        {
            (*this).erase(*pos);
        }

        return (*this);

    }  // End of method 'typeValueMap::removeTypeID()'



      /* Return the data value (double) corresponding to provided type.
       *
       * @param type       Type of value to be looked for.
       */
    double typeValueMap::getValue(const TypeID& type) const
        throw(TypeIDNotFound)
    {

        typeValueMap::const_iterator itObs( (*this).find(type) );
        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(TypeIDNotFound("TypeID not found in map"));
        }

    }  // End of method 'typeValueMap::getValue()'



      // Return a reference to the data value (double) with
      // corresponding type.
      // @param type Type of value to be looked for.
    double& typeValueMap::operator()(const TypeID& type)
        throw(TypeIDNotFound)
    {

        typeValueMap::iterator itObs ( (*this).find(type) );

        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(TypeIDNotFound("TypeID not found in map"));
        }

    }  // End of method 'typeValueMap::operator()'



      ////// satValueMap //////


      // Return a SatIDSet with all the satellites present in this object.
    SatIDSet satValueMap::getSatID() const
    {

        SatIDSet satSet;

        for( satValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            satSet.insert( (*pos).first );
        }

        return satSet;

    }  // End of method 'satValueMap::getSatID()'



      // Return a Vector with all the satellites present in this object.
    Vector<SatID> satValueMap::getVectorOfSatID() const
    {

        std::vector<SatID> temp;

        for( satValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            temp.push_back( (*pos).first );
        }

        Vector<SatID> result;
        result = temp;

        return result;

    }  // End of method 'satValueMap::getVectorOfSatID()'



      // Return a satValueMap with only this satellite.
      // @param satellite Satellite to be extracted.
    satValueMap satValueMap::extractSatID(const SatID& satellite) const
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return extractSatID(satSet);

    }  // End of method 'satValueMap::extractSatID()'



      // Return a satValueMap with only one satellite, identified by
      // the given parameters.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    satValueMap satValueMap::extractSatID( const int& p,
                                           const SatID::SatelliteSystem& s ) const
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return extractSatID(tempSatellite);

    }  // End of method 'satValueMap::extractSatID()'



      // Return a satValueMap with only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to
      //               be extracted.
    satValueMap satValueMap::extractSatID(const SatIDSet& satSet) const
    {

        satValueMap svMap;

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            satValueMap::const_iterator itObs( (*this).find(*pos) );

            if( itObs != (*this).end() )
            {
                svMap[ (*itObs).first ] = (*itObs).second;
            }
        }

        return svMap;

    }  // End of method 'satValueMap::extractSatID()'



      // Modifies this object, keeping only this satellite.
      // @param satellite Satellite to be kept.
    satValueMap& satValueMap::keepOnlySatID(const SatID& satellite)
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return keepOnlySatID(satSet);

    }  // End of method 'satValueMap::keepOnlySatID()'



      // Modifies this object, keeping only this satellite.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    satValueMap& satValueMap::keepOnlySatID( const int& p,
                                             const SatID::SatelliteSystem& s )
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return keepOnlySatID(tempSatellite);

    }  // End of method 'satValueMap::keepOnlySatID()'



      // Modifies this object, keeping only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to be kept.
    satValueMap& satValueMap::keepOnlySatID(const SatIDSet& satSet)
    {

        satValueMap svMap = extractSatID(satSet);
        (*this) = svMap;

        return (*this);

    }  // End of method 'satValueMap::keepOnlySatID()'



      // Modifies this object, removing these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to
      //               be removed.
    satValueMap& satValueMap::removeSatID(const SatIDSet& satSet)
    {

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            (*this).erase(*pos);
        }

        return (*this);

    }  // End of method 'satValueMap::removeSatID()'



      /* Return the data value (double) corresponding to provided SatID.
       *
       * @param satellite     Satellite to be looked for.
       */
    double satValueMap::getValue(const SatID& satellite) const
        throw(SatIDNotFound)
    {

        satValueMap::const_iterator itObs( (*this).find(satellite) );
        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(SatIDNotFound("SatID not found in map"));
        }

    }  // End of method 'satValueMap::getValue()'



      // Return a reference to the data value (double) with
      // corresponding SatID.
      // @param satellite Satellite to be looked for.
    double& satValueMap::operator()(const SatID& satellite)
        throw(SatIDNotFound)
    {

        satValueMap::iterator itObs( (*this).find(satellite) );

        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(SatIDNotFound("SatID not found in map"));
        }

    }  // End of method 'satValueMap::operator()'



      ////// sourceValueMap //////


      // Return a SourceIDSet with all the sources present in this object.
    SourceIDSet sourceValueMap::getSourceID() const
    {

        SourceIDSet sourceSet;

        for( sourceValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            sourceSet.insert( (*pos).first );
        }

        return sourceSet;

    }  // End of method 'sourceValueMap::getSourceID()'



      // Return a Vector with all the sources present in this object.
    Vector<SourceID> sourceValueMap::getVectorOfSourceID() const
    {

        std::vector<SourceID> temp;

        for( sourceValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            temp.push_back( (*pos).first );
        }

        Vector<SourceID> result;
        result = temp;

        return result;

    }  // End of method 'sourceValueMap::getVectorOfSourceID()'



      // Return a sourceValueMap with only this source.
      // @param source Source to be extracted.
    sourceValueMap sourceValueMap::extractSourceID(const SourceID& source) const
    {

        SourceIDSet sourceSet;
        sourceSet.insert(source);

        return extractSourceID(sourceSet);

    }  // End of method 'sourceValueMap::extractSourceID()'



      // Return a sourceValueMap with only these sources.
      // @param sourceSet Set (SourceIDSet) containing the sources to
      //               be extracted.
    sourceValueMap sourceValueMap::extractSourceID(const SourceIDSet& sourceSet) const
    {

        sourceValueMap svMap;

        for( SourceIDSet::const_iterator pos = sourceSet.begin();
             pos != sourceSet.end();
             ++pos )
        {
            sourceValueMap::const_iterator itObs( (*this).find(*pos) );

            if( itObs != (*this).end() )
            {
                svMap[ (*itObs).first ] = (*itObs).second;
            }
        }

        return svMap;

    }  // End of method 'sourceValueMap::extractSourceID()'



      // Modifies this object, keeping only this source.
      // @param source Source to be kept.
    sourceValueMap& sourceValueMap::keepOnlySourceID(const SourceID& source)
    {

        SourceIDSet sourceSet;
        sourceSet.insert(source);

        return keepOnlySourceID(sourceSet);

    }  // End of method 'sourceValueMap::keepOnlySourceID()'



      // Modifies this object, keeping only these sources.
      // @param sourceSet Set (SourceIDSet) containing the sources to be kept.
    sourceValueMap& sourceValueMap::keepOnlySourceID(const SourceIDSet& sourceSet)
    {

        sourceValueMap svMap = extractSourceID(sourceSet);
        (*this) = svMap;

        return (*this);

    }  // End of method 'sourceValueMap::keepOnlySourceID()'



      // Modifies this object, removing these sources.
      // @param sourceSet Set (SourceIDSet) containing the sources to
      //               be removed.
    sourceValueMap& sourceValueMap::removeSourceID(const SourceIDSet& sourceSet)
    {

        for( SourceIDSet::const_iterator pos = sourceSet.begin();
             pos != sourceSet.end();
             ++pos )
        {
            (*this).erase(*pos);
        }

        return (*this);

    }  // End of method 'sourceValueMap::removeSourceID()'



      /* Return the data value (double) corresponding to provided SourceID.
       *
       * @param source     Source to be looked for.
       */
    double sourceValueMap::getValue(const SourceID& source) const
        throw(SourceIDNotFound)
    {

        sourceValueMap::const_iterator itObs( (*this).find(source) );
        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(SourceIDNotFound("SourceID not found in map"));
        }

    }  // End of method 'sourceValueMap::getValue()'



      // Return a reference to the data value (double) with
      // corresponding SourceID.
      // @param source Source to be looked for.
    double& sourceValueMap::operator()(const SourceID& source)
        throw(SourceIDNotFound)
    {

        sourceValueMap::iterator itObs( (*this).find(source) );

        if ( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(SourceIDNotFound("SourceID not found in map"));
        }

    }  // End of method 'sourceValueMap::operator()'



      ////// satTypeValueMap //////


      /* Return the total number of data elements in the map.
       * This method DOES NOT suppose that all the satellites have
       * the same number of type values.
       */
    size_t satTypeValueMap::numElements() const
    {

        size_t numEle(0);

        for( satTypeValueMap::const_iterator it = (*this).begin();
             it != (*this).end();
             ++it )
        {
            numEle = numEle + (*it).second.size();
        }

        return numEle;

    }  // End of method 'satTypeValueMap::numElements()'



      // Return a SatIDSet with all the satellites present in this object.
    SatIDSet satTypeValueMap::getSatID() const
    {

        SatIDSet satSet;

        for( satTypeValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            satSet.insert( (*pos).first );
        }

        return satSet;

    }  // End of method 'satTypeValueMap::getSatID()'



      // Return a Vector with all the satellites present in this object.
    Vector<SatID> satTypeValueMap::getVectorOfSatID() const
    {

        std::vector<SatID> temp;

        for( satTypeValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {
            temp.push_back( (*pos).first );
        }

        Vector<SatID> result;
        result = temp;

        return result;

    }  // End of method 'satTypeValueMap::getVectorOfSatID()'



      // Return a TypeIDSet with all the data types present in
      // this object. This does not imply that all satellites have these types.
    TypeIDSet satTypeValueMap::getTypeID() const
    {

        TypeIDSet typeSet;

        for( satTypeValueMap::const_iterator pos = (*this).begin();
             pos != (*this).end();
             ++pos )
        {

            for( typeValueMap::const_iterator it = (*pos).second.begin();
                 it != (*pos).second.end();
                 ++it )
            {
                typeSet.insert( (*it).first );
            }

        }

        return typeSet;

    }  // End of method 'satTypeValueMap::getTypeID()'



      // Return a satTypeValueMap with only this satellite.
      // @param satellite Satellite to be extracted.
    satTypeValueMap satTypeValueMap::extractSatID(const SatID& satellite) const
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return extractSatID(satSet);

    }  // End of method 'satTypeValueMap::extractSatID()'



      // Return a satTypeValueMap with only one satellite, identified
      // by the given parameters.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    satTypeValueMap satTypeValueMap::extractSatID( const int& p,
                                                   const SatID::SatelliteSystem& s ) const
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return extractSatID(tempSatellite);

    }  // End of method 'satTypeValueMap::extractSatID()'



      // Return a satTypeValueMap with only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to
      //               be extracted.
    satTypeValueMap satTypeValueMap::extractSatID(const SatIDSet& satSet) const
    {

        satTypeValueMap stvMap;

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            satTypeValueMap::const_iterator itObs( (*this).find(*pos) );
            if( itObs != (*this).end() )
            {
                stvMap[ (*itObs).first ] = (*itObs).second;
            }
        }

        return stvMap;

    }  // End of method 'satTypeValueMap::extractSatID()'



      // Modifies this object, keeping only this satellite.
      // @param satellite Satellite to be kept.
    satTypeValueMap& satTypeValueMap::keepOnlySatID(const SatID& satellite)
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return keepOnlySatID(satSet);

    }  // End of method 'satTypeValueMap::keepOnlySatID()'



      // Modifies this object, keeping only this satellite.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    satTypeValueMap& satTypeValueMap::keepOnlySatID( const int& p,
                                                     const SatID::SatelliteSystem& s )
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return keepOnlySatID(tempSatellite);

    }  // End of method 'satTypeValueMap::keepOnlySatID()'



      // Modifies this object, keeping only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to be kept.
    satTypeValueMap& satTypeValueMap::keepOnlySatID(const SatIDSet& satSet)
    {

        satTypeValueMap stvMap( extractSatID(satSet) );
        (*this) = stvMap;

        return (*this);

    }  // End of method 'satTypeValueMap::keepOnlySatID()'



      // Return a satTypeValueMap with only this type of value.
      // @param type Type of value to be extracted.
    satTypeValueMap satTypeValueMap::extractTypeID(const TypeID& type) const
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return extractTypeID(typeSet);

    }  // End of method 'satTypeValueMap::extractTypeID()'



      // Return a satTypeValueMap with only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be extracted.
    satTypeValueMap satTypeValueMap::extractTypeID(const TypeIDSet& typeSet) const
    {

        satTypeValueMap theMap;

        for( satTypeValueMap::const_iterator it = (*this).begin();
             it != (*this).end();
             ++it )
        {

            typeValueMap tvMap( (*it).second.extractTypeID(typeSet) );
            if( tvMap.size() > 0 )
            {
                theMap[(*it).first] = tvMap;
            }

        }

        return theMap;

    }  // End of method 'satTypeValueMap::extractTypeID()'



      // Modifies this object, keeping only this type of data.
      // @param type Type of value to be kept.
    satTypeValueMap& satTypeValueMap::keepOnlyTypeID(const TypeID& type)
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return keepOnlyTypeID(typeSet);

    }  // End of method 'satTypeValueMap::keepOnlyTypeID()'



      // Modifies this object, keeping only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    satTypeValueMap& satTypeValueMap::keepOnlyTypeID(const TypeIDSet& typeSet)
    {

        satTypeValueMap stvMap( extractTypeID(typeSet) );
        (*this) = stvMap;

        return (*this);

    }  // End of method 'satTypeValueMap::keepOnlyTypeID()'



      // Modifies this object, removing these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be removed.
    satTypeValueMap& satTypeValueMap::removeSatID(const SatIDSet& satSet)
    {

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            (*this).erase(*pos);
        }

        return (*this);

    }  // End of method 'satTypeValueMap::removeSatID()'



      // Modifies this object, removing this type of data.
      // @param type Type of value to be removed.
    satTypeValueMap& satTypeValueMap::removeTypeID(const TypeID& type)
    {

        for( satTypeValueMap::iterator it = (*this).begin();
             it != (*this).end();
             ++it )
        {
            (*it).second.removeTypeID(type);
        }

        return (*this);

    }  // End of method 'satTypeValueMap::removeTypeID()'



      // Modifies this object, removing these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    satTypeValueMap& satTypeValueMap::removeTypeID(const TypeIDSet& typeSet)
    {

        for( TypeIDSet::const_iterator pos = typeSet.begin();
             pos != typeSet.end();
             ++pos )
        {
            removeTypeID(*pos);
        }

        return (*this);

    }  // End of method 'satTypeValueMap::removeTypeID()'



      // Return a GPSTk::Vector containing the data values with this type.
      // @param type Type of value to be returned.
      // This method returns zero if a given satellite does not have this type.
    Vector<double> satTypeValueMap::getVectorOfTypeID(const TypeID& type) const
    {

         // Let's declare a STL vector
        std::vector<double> temp;

        for( satTypeValueMap::const_iterator it = (*this).begin();
             it != (*this).end();
             ++it )
        {

            typeValueMap::const_iterator itObs( (*it).second.find(type) );
            if ( itObs != (*it).second.end() )
            {
                temp.push_back( (*itObs).second );
            }
            else
            {
                temp.push_back( 0.0 );
            }

        }

         // Let's declare a GPSTk Vector
        Vector<double> result;

         // Transform STL vector into GPSTk Vector
        result = temp;

        return result;

    }  // End of method 'satTypeValueMap::getVectorOfTypeID()'



      // Return a GPSTk::Matrix containing the data values in this set.
      // @param typeSet  TypeIDSet of values to be returned.
    Matrix<double> satTypeValueMap::getMatrixOfTypes(const TypeIDSet& typeSet) const
    {

         // First, let's create a Matrix<double> of the proper size
        Matrix<double> tempMat( (*this).numSats(), typeSet.size(), 0.0 );

        size_t numRow(0), numCol(0);

        for( satTypeValueMap::const_iterator it = (*this).begin();
             it != (*this).end();
             ++it )
        {
            numCol=0;

            for( TypeIDSet::const_iterator pos = typeSet.begin();
                 pos != typeSet.end();
                 ++pos )
            {

                typeValueMap::const_iterator itObs( (*it).second.find(*pos) );
                if( itObs != (*it).second.end() )
                {
                    tempMat(numRow, numCol) = (*itObs).second;
                }

                ++numCol;
            }

            ++numRow;

        }

        return tempMat;

    }  // End of method 'satTypeValueMap::getMatrixOfTypes()'



      /* Modifies this object, adding one vector of data with this type,
       * one value per satellite.
       *
       * If type already exists, data is overwritten. If the number of
       * values does not match with the number of satellites, a
       * NumberOfSatsMismatch exception is thrown.
       *
       * Given that dataVector does not store information about the
       * satellites the values correspond to, the user is held responsible
       * for having the data values stored in dataVector in the proper
       * order regarding the SatIDs in this object.
       *
       * @param type          Type of data to be added.
       * @param dataVector    GPSTk Vector containing the data to be added.
       */
    satTypeValueMap& satTypeValueMap::insertTypeIDVector( const TypeID& type,
                                                          const Vector<double> dataVector )
        throw(NumberOfSatsMismatch)
    {

        if( dataVector.size() == (*this).numSats() )
        {
            size_t pos = 0;

            for( satTypeValueMap::iterator it = (*this).begin();
                 it != (*this).end();
                 ++it )
            {
                (*it).second[type] = dataVector[pos];
                ++pos;
            }

            return (*this);

        }
        else
        {
            GPSTK_THROW( NumberOfSatsMismatch(" Number of data values in vector \
and number of satellites do not match") );
        }

    }  // End of method 'satTypeValueMap::insertTypeIDVector()'



      /* Modifies this object, adding a matrix of data, one vector
       * per satellite.
       *
       * If types already exists, data is overwritten. If the number of
       * rows in matrix does not match with the number of satellites, a
       * NumberOfSatsMismatch exception is thrown. If the number of columns
       * in matrix does not match with the number of types in typeSet, a
       * NumberOfTypesMismatch exception is thrown.
       *
       * Given that dataMatrix does not store information about the
       * satellites and types the values correspond to, the user is held
       * responsible for having those data values stored in dataMatrix in
       * the proper order regarding the SatIDs in this object and the
       * provided typeSet.
       *
       * @param typeSet       Set (TypeIDSet) containing the types of data
       *                      to be added.
       * @param dataMatrix    GPSTk Matrix containing the data to be added.
       */
    satTypeValueMap& satTypeValueMap::insertMatrix( const TypeIDSet& typeSet,
                                                    const Matrix<double> dataMatrix )
        throw(NumberOfSatsMismatch, NumberOfTypesMismatch)
    {

        if( dataMatrix.rows() != (*this).numSats() )
        {
            GPSTK_THROW( NumberOfSatsMismatch("Number of rows in matrix and \
number of satellites do not match") );
        }

        if( dataMatrix.cols() == typeSet.size() )
        {

            size_t pos(0);

            for( satTypeValueMap::iterator it = (*this).begin();
                 it != (*this).end();
                 ++it )
            {

                size_t idx(0);

                for( TypeIDSet::const_iterator itSet = typeSet.begin();
                     itSet != typeSet.end();
                     ++itSet )
                {
                    (*it).second[(*itSet)] = dataMatrix(pos,idx);
                    ++idx;
                }

                ++pos;

            }

            return (*this);

        }
        else
        {
            GPSTK_THROW( NumberOfTypesMismatch("Number of data values per row \
in matrix and number of types do not match") );
        }

    }  // End of method 'satTypeValueMap::insertMatrix()'



      /* Return the data value (double) corresponding to provided SatID
       * and TypeID.
       *
       * @param satellite     Satellite to be looked for.
       * @param type          Type to be looked for.
       */
    double satTypeValueMap::getValue( const SatID& satellite,
                                      const TypeID& type ) const
        throw( SatIDNotFound, TypeIDNotFound )
    {

        satTypeValueMap::const_iterator itObs( (*this).find(satellite) );
        if( itObs != (*this).end() )
        {
            return (*itObs).second.getValue( type );
        }
        else
        {
            GPSTK_THROW(SatIDNotFound("SatID not found in map"));
        }

    }  // End of method 'satTypeValueMap::getValue()'



      // Return a reference to the typeValueMap with corresponding SatID.
      // @param type Type of value to be looked for.
    typeValueMap& satTypeValueMap::operator()(const SatID& satellite)
        throw(SatIDNotFound)
    {

        satTypeValueMap::iterator itObs( (*this).find(satellite) );
        if( itObs != (*this).end() )
        {
            return (*itObs).second;
        }
        else
        {
            GPSTK_THROW(SatIDNotFound("SatID not found in map"));
        }

    }  // End of method 'satTypeValueMap::operator()'



      ////// gnssSatValue //////


      // Return a gnssSatValue with only this satellite.
      // @param satellite Satellite to be extracted.
    gnssSatValue gnssSatValue::extractSatID(const SatID& satellite) const
    {

        gnssSatValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satellite);

        return result;

    }  // End of method 'gnssSatValue::extractSatID()'



      // Return a gnssSatValue with only one satellite, identified by
      // the given parameters.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssSatValue gnssSatValue::extractSatID( const int& p,
                                             const SatID::SatelliteSystem& s ) const
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return extractSatID(tempSatellite);

    }  // End of method 'gnssSatValue::extractSatID()'



      // Return a gnssSatValue with only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be extracted.
    gnssSatValue gnssSatValue::extractSatID(const SatIDSet& satSet) const
    {

        gnssSatValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satSet);

        return result;

    }  // End of method 'gnssSatValue::extractSatID()'



      // Modifies this object, keeping only this satellite.
      // @param satellite Satellite to be kept.
    gnssSatValue& gnssSatValue::keepOnlySatID(const SatID& satellite)
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return keepOnlySatID(satSet);

    }  // End of method 'gnssSatValue::keepOnlySatID()'



      // Modifies this object, keeping only this satellite.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssSatValue& gnssSatValue::keepOnlySatID( const int& p,
                                               const SatID::SatelliteSystem& s )
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return keepOnlySatID(tempSatellite);

    }  // End of method 'gnssSatValue::keepOnlySatID()'



      // Modifies this object, keeping only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to be kept.
    gnssSatValue& gnssSatValue::keepOnlySatID(const SatIDSet& satSet)
    {

        satValueMap svMap ( (*this).body.extractSatID(satSet) );
        (*this).body = svMap;

        return (*this);

    }  // End of method 'gnssSatValue::keepOnlySatID()'



      // Modifies this object, removing these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be removed.
    gnssSatValue& gnssSatValue::removeSatID(const SatIDSet& satSet)
    {

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            (*this).body.erase(*pos);
        }

        return (*this);

    }  // End of method 'gnssSatValue::removeSatID()'




      ////// gnssTypeValue //////


      // Return a gnssTypeValue with only this type of data.
      // @param type Type of value to be extracted.
    gnssTypeValue gnssTypeValue::extractTypeID(const TypeID& type) const
    {

        gnssTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(type);

        return result;

    }  // End of method 'gnssTypeValue::extractTypeID()'



      // Return a gnssTypeValue with only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be extracted.
    gnssTypeValue gnssTypeValue::extractTypeID(const TypeIDSet& typeSet) const
    {

        gnssTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(typeSet);

        return result;

    }  // End of method 'gnssTypeValue::extractTypeID()'



      // Modifies this object, keeping only this type of data.
      // @param type Type of value to be kept.
    gnssTypeValue& gnssTypeValue::keepOnlyTypeID(const TypeID& type)
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return keepOnlyTypeID(typeSet);

    }  // End of method 'gnssTypeValue::keepOnlyTypeID()'



      // Modifies this object, keeping only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    gnssTypeValue& gnssTypeValue::keepOnlyTypeID(const TypeIDSet& typeSet)
    {

        typeValueMap tvMap( (*this).body.extractTypeID(typeSet) );
        (*this).body = tvMap;

        return (*this);

    }  // End of method 'gnssTypeValue::keepOnlyTypeID()'



      // Modifies this object, removing these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    gnssTypeValue& gnssTypeValue::removeTypeID(const TypeIDSet& typeSet)
    {

        for( TypeIDSet::const_iterator pos = typeSet.begin();
             pos != typeSet.end();
             ++pos )
        {
            (*this).body.erase(*pos);
        }

        return (*this);

    }  // End of method 'gnssTypeValue::removeTypeID()'




      ////// gnssSatTypeValue //////


      // Return a gnssSatTypeValue with only this satellite.
      // @param satellite Satellite to be extracted.
    gnssSatTypeValue gnssSatTypeValue::extractSatID(const SatID& satellite) const
    {

        gnssSatTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satellite);

        return result;

    }  // End of method 'gnssSatTypeValue::extractSatID()'



      // Return a gnssSatTypeValue with only one satellite, identified
      // by the given parameters.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssSatTypeValue gnssSatTypeValue::extractSatID( const int& p,
                                                     const SatID::SatelliteSystem& s ) const
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return extractSatID(tempSatellite);

    }  // End of method 'gnssSatTypeValue::extractSatID()'



      // Return a gnssSatTypeValue with only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be extracted.
    gnssSatTypeValue gnssSatTypeValue::extractSatID(const SatIDSet& satSet) const
    {

        gnssSatTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satSet);

        return result;

    }  // End of method 'gnssSatTypeValue::extractSatID()'



      // Modifies this object, keeping only this satellite.
      // @param satellite Satellite to be kept.
    gnssSatTypeValue& gnssSatTypeValue::keepOnlySatID(const SatID& satellite)
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return keepOnlySatID(satSet);

    }  // End of method 'gnssSatTypeValue::keepOnlySatID()'



      // Modifies this object, keeping only this satellite.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssSatTypeValue& gnssSatTypeValue::keepOnlySatID( const int& p,
                                                       const SatID::SatelliteSystem& s )
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return keepOnlySatID(tempSatellite);

    }  // End of method 'gnssSatTypeValue::keepOnlySatID()'



      // Modifies this object, keeping only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to be kept.
    gnssSatTypeValue& gnssSatTypeValue::keepOnlySatID(const SatIDSet& satSet)
    {

        satTypeValueMap stvMap( (*this).body.extractSatID(satSet) );
        (*this).body = stvMap;

        return (*this);

    }  // End of method 'gnssSatTypeValue::keepOnlySatID()'



      // Return a gnssSatTypeValue with only this type of data.
      // @param type Type of value to be extracted.
    gnssSatTypeValue gnssSatTypeValue::extractTypeID(const TypeID& type) const
    {

        gnssSatTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(type);

        return result;

    }  // End of method 'gnssSatTypeValue::extractTypeID()'



      // Return a gnssSatTypeValue with only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be extracted.
    gnssSatTypeValue gnssSatTypeValue::extractTypeID(const TypeIDSet& typeSet) const
    {

        gnssSatTypeValue result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(typeSet);

        return result;

    }  // End of method 'gnssSatTypeValue::extractTypeID()'



      // Modifies this object, keeping only this type of data.
      // @param type Type of value to be kept.
    gnssSatTypeValue& gnssSatTypeValue::keepOnlyTypeID(const TypeID& type)
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return keepOnlyTypeID(typeSet);

    }  // End of method 'gnssSatTypeValue::keepOnlyTypeID()'



      // Modifies this object, keeping only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    gnssSatTypeValue& gnssSatTypeValue::keepOnlyTypeID(const TypeIDSet& typeSet)
    {

        satTypeValueMap stvMap( (*this).body.extractTypeID(typeSet) );
        (*this).body = stvMap;

        return (*this);

    }  // End of method 'gnssSatTypeValue::keepOnlyTypeID()'



      // Modifies this object, removing these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be removed.
    gnssSatTypeValue& gnssSatTypeValue::removeSatID(const SatIDSet& satSet)
    {

        for( SatIDSet::const_iterator pos = satSet.begin();
             pos != satSet.end();
             ++pos )
        {
            (*this).body.erase(*pos);
        }

        return (*this);

    }  // End of method 'gnssSatTypeValue::removeSatID()'



      // Modifies this object, removing these types of data
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    gnssSatTypeValue& gnssSatTypeValue::removeTypeID(const TypeIDSet& typeSet)
    {

        for( TypeIDSet::const_iterator pos = typeSet.begin();
             pos != typeSet.end();
             ++pos )
        {
            (*this).body.removeTypeID(*pos);
        }

        return (*this);

    }  // End of method 'gnssSatTypeValue::removeTypeID()'




      ////// gnssRinex //////


      // Return a gnssRinex with only this satellite.
      // @param satellite Satellite to be extracted.
    gnssRinex gnssRinex::extractSatID(const SatID& satellite) const
    {

        gnssRinex result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satellite);

        return result;

    }  // End of method 'gnssRinex::extractSatID()'



      // Return a gnssRinex with only one satellite, identified by
      // the given parameters.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssRinex gnssRinex::extractSatID( const int& p,
                                       const SatID::SatelliteSystem& s ) const
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return extractSatID(tempSatellite);

    }  // End of method 'gnssRinex::extractSatID()'



      // Return a gnssRinex with only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites
      //               to be extracted.
    gnssRinex gnssRinex::extractSatID(const SatIDSet& satSet) const
    {

        gnssRinex result;
        result.header = (*this).header;
        result.body = (*this).body.extractSatID(satSet);

        return result;

    }  // End of method 'gnssRinex::extractSatID()'



      // Modifies this object, keeping only this satellite.
      // @param satellite Satellite to be kept.
    gnssRinex& gnssRinex::keepOnlySatID(const SatID& satellite)
    {

        SatIDSet satSet;
        satSet.insert(satellite);

        return keepOnlySatID(satSet);

    }  // End of method 'gnssRinex::keepOnlySatID()'



      // Modifies this object, keeping only this satellite.
      // @param p Satellite PRN number.
      // @param p System the satellite belongs to.
    gnssRinex& gnssRinex::keepOnlySatID( const int& p,
                                         const SatID::SatelliteSystem& s )
    {

        SatID tempSatellite(p, s);  // We build a temporary SatID object

        return keepOnlySatID(tempSatellite);

    }  // End of method 'gnssRinex::keepOnlySatID()'



      // Modifies this object, keeping only these satellites.
      // @param satSet Set (SatIDSet) containing the satellites to be kept.
    gnssRinex& gnssRinex::keepOnlySatID(const SatIDSet& satSet)
    {

        satTypeValueMap stvMap( (*this).body.extractSatID(satSet) );
        (*this).body = stvMap;

        return (*this);

    }  // End of method 'gnssRinex::keepOnlySatID()'



      // Return a gnssRinex with only this type of data.
      // @param type Type of value to be extracted.
    gnssRinex gnssRinex::extractTypeID(const TypeID& type) const
    {

        gnssRinex result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(type);

        return result;

    }  // End of method 'gnssRinex::extractTypeID()'



      // Return a gnssRinex with only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be extracted.
    gnssRinex gnssRinex::extractTypeID(const TypeIDSet& typeSet) const
    {

        gnssRinex result;
        result.header = (*this).header;
        result.body = (*this).body.extractTypeID(typeSet);

        return result;

    }  // End of method 'gnssRinex::extractTypeID()'



      // Modifies this object, keeping only this type of data.
      // @param type Type of value to be kept.
    gnssRinex& gnssRinex::keepOnlyTypeID(const TypeID& type)
    {

        TypeIDSet typeSet;
        typeSet.insert(type);

        return keepOnlyTypeID(typeSet);

    }  // End of method 'gnssRinex::keepOnlyTypeID()'



      // Modifies this object, keeping only these types of data.
      // @param typeSet Set (TypeIDSet) containing the types of data
      //                to be kept.
    gnssRinex& gnssRinex::keepOnlyTypeID(const TypeIDSet& typeSet)
    {

        satTypeValueMap stvMap( (*this).body.extractTypeID(typeSet) );
        (*this).body = stvMap;

        return (*this);

    }  // End of method 'gnssRinex::keepOnlyTypeID()'



      // Return a gnssRinex with only these types of data.
      // @param satSys Satellite System value to be kept.
    gnssRinex& gnssRinex::keepOnlySatSystem(const SatID::SatelliteSystem satSys)
    {
        satTypeValueMap stvMap( (*this).body );

        SatIDSet satRejectedSet;
        for(satTypeValueMap::iterator it = stvMap.begin();
            it != stvMap.end();
            ++it)
        {
            if( it->first.system != satSys ) satRejectedSet.insert( it->first );
        }
        stvMap.removeSatID(satRejectedSet);


        (*this).body = stvMap;

        return (*this);

    }  // End of method 'gnssRinex::keepOnlySatSystem()'



      //// Some other handy data structures



      /* Return the data value (double) corresponding to provided SourceID,
       * SatID and TypeID.
       *
       * @param source        Source to be looked for.
       * @param satellite     Satellite to be looked for.
       * @param type          Type to be looked for.
       */
    double sourceDataMap::getValue( const SourceID& source,
                                    const SatID& satellite,
                                    const TypeID& type ) const
        throw( SourceIDNotFound, SatIDNotFound, TypeIDNotFound )
    {

         // Look for the SourceID
        sourceDataMap::const_iterator itObs( (*this).find(source) );
        if( itObs != (*this).end() )
        {
            return (*itObs).second.getValue( satellite, type );
        }
        else
        {
            GPSTK_THROW(SourceIDNotFound("SourceID not found in map"));
        }

    }  // End of method 'sourceDataMap::getValue()'



      /* Get a set with all the SourceID's in this data structure.
       *
       * @warning If current 'sourceDataMap' is big, this could be a very
       * costly operation.
       */
    SourceIDSet sourceDataMap::getSourceIDSet( void ) const
    {

         // SourceID set to be returned
        SourceIDSet toReturn;

         // Then, iterate through corresponding 'sourceDataMap'
        for(sourceDataMap::const_iterator sdmIter = (*this).begin();
            sdmIter != (*this).end();
            ++sdmIter)
        {

            // Add current SourceID to 'toReturn'
            toReturn.insert( (*sdmIter).first );

        }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        return toReturn;

    }  // End of method 'sourceDataMap::getSourceIDSet()'



      /* Get a set with all the SatID's in this data structure.
       *
       * @warning If current 'sourceDataMap' is big, this could be a very
       * costly operation.
       */
    SatIDSet sourceDataMap::getSatIDSet( void ) const
    {

         // SatID set to be returned
        SatIDSet toReturn;

         // Then, iterate through corresponding 'sourceDataMap'
        for(sourceDataMap::const_iterator sdmIter = (*this).begin();
            sdmIter != (*this).end();
            ++sdmIter)
        {

            // Finally, iterate through corresponding 'satTypeValueMap'
            for(satTypeValueMap::const_iterator stvmIter =
                                                   (*sdmIter).second.begin();
                stvmIter != (*sdmIter).second.end();
                stvmIter++)
            {

                // Add current SatID to 'toReturn'
                toReturn.insert( (*stvmIter).first );

            }  // End of 'for( satTypeValueMap::const_iterator stvmIter = ...'

        }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        return toReturn;

    }  // End of method 'sourceDataMap::getSatIDSet()'



      /* Adds 'gnssSatTypeValue' object data to this structure.
       *
       * @param gds     gnssSatTypeValue object containing data to be added.
       */
    gnssDataMap& gnssDataMap::addGnssSatTypeValue( const gnssSatTypeValue& gds )
    {

         // Declare an object for temporal storage of data
        sourceDataMap sdMap;

         // Fill with data
        sdMap[ gds.header.source ] = gds.body;

         // Introduce data into this GDS
        (*this).insert( pair<const CommonTime, sourceDataMap>( gds.header.epoch, sdMap ) );

         // Return curren GDS.
        return (*this);

    }  // End of method 'gnssDataMap::addGnssSatTypeValue()'



      /* Adds 'gnssRinex' object data to this structure.
       *
       * @param gds     gnssRinex object containing data to be added.
       */
    gnssDataMap& gnssDataMap::addGnssRinex( const gnssRinex& gds )
    {

         // Declare an object for temporal storage of data
        sourceDataMap sdMap;

         // Fill with data
        sdMap[ gds.header.source ] = gds.body;

         // Introduce data into this GDS
        (*this).insert( pair<const CommonTime, sourceDataMap>( gds.header.epoch, sdMap ) );

         // Return current GDS.
        return (*this);

    }  // End of method 'gnssDataMap::addGnssRinex()'



      /* Return a 'gnssRinex' object corresponding to given SourceID.
       *
       * @param source     SourceID object.
       *
       * @warning Returned data will correspond to first matching SourceID,
       * if it exists.
       */
    gnssRinex gnssDataMap::getGnssRinex( const SourceID& source ) const
    {

         // Get first data set
        gnssDataMap gdMap( (*this).frontEpoch() );

         // Declare a gnssRinex object to be returned
        gnssRinex toReturn;

         // We'll need a flag
        bool found(false);

         // Look into the data structure
        for(gnssDataMap::const_iterator it = gdMap.begin();
            it != gdMap.end() && !found;
            ++it)
        {
            // Always return the epoch/epochFlag
            toReturn.header.epoch = (*it).first;
            toReturn.header.epochFlag = 0;

            // Look for SourceID
            sourceDataMap::const_iterator iter( (*it).second.find( source ) );
            if( iter != (*it).second.end() )
            {
                toReturn.body = (*iter).second;
                toReturn.header.source = (*iter).first;
                found = true;
            }

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

         // Return GDS
        return toReturn;

    }  // End of method 'gnssDataMap::getGnssRinex()'


      /** Return a 'gnssRinex' object corresponding to given SatID.
       *
       * @param sate       SatID object.
       *
       * @warning The first gnssRinex which contains the given SatID
       * is returned.
       */
    gnssRinex gnssDataMap::getGnssRinex( const SatID& sat ) const
    {

         // Get first data set
        gnssDataMap gdMap( (*this).frontEpoch() );

         // Declare a gnssRinex object to be returned
        gnssRinex toReturn;

         // We'll need a flag
        bool found(false);

         // Look into the data structure
        for(gnssDataMap::const_iterator it = gdMap.begin();
            it != gdMap.end() && !found;
            ++it)
        {

            // Then, iterate through corresponding 'sourceDataMap'
            for(sourceDataMap::const_iterator sdmIter = (*it).second.begin();
                sdmIter != (*it).second.end() && !found;
                ++sdmIter)
            {
                // Look for the satellite
                satTypeValueMap::const_iterator stmIter( (*sdmIter).second.find(sat) );

                // Found
                if( stmIter != (*sdmIter).second.end() )
                {
                    // Get the gnssRinex related to this source
                    toReturn = getGnssRinex( (*sdmIter).first );

                    // Work is done, let's get out
                    found = true;

                }  // End of 'if( stmIter != ... )'

            }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        // Return GDS
        return toReturn;

    }  // End of method 'gnssDataMap::getGnssRinex()'


      /* Adds 'gnssDataMap' object data to this structure.
       *
       * @param gds     gnssDataMap object containing data to be added.
       */
    gnssDataMap& gnssDataMap::addGnssDataMap( const gnssDataMap& gds )
    {
         // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = gds.begin();
            it != gds.end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue stv;
                stv.header.epoch = time;
                stv.header.source = itsrc->first;
                stv.body = itsrc->second;

                addGnssSatTypeValue(stv);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return (*this);
    }


      // Return a copy of the first element in the map.
    gnssDataMap gnssDataMap::front( void ) const
    {

         // Object to be returned
        gnssDataMap toReturn;

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get first element
            gnssDataMap::const_iterator firstIt( (*this).begin() );

            // Insert first element in 'toReturn' object
            toReturn.insert(*firstIt);

        }

        return toReturn;

    }  // End of method 'gnssDataMap::front()'



      // Removes the first element in the map.
    void gnssDataMap::pop_front( void )
    {

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Define an iterator pointing to the beginning of this map
            gnssDataMap::iterator pos( (*this).begin() );

            // It is advisable to avoid sawing off the branch we are sitting on
            (*this).erase( pos++ );

        }

        return;

    }  // End of method 'gnssDataMap::pop_front()'


      // Return the data corresponding to the first epoch in the map,
      // taking into account the 'tolerance' value.
    gnssDataMap gnssDataMap::frontEpoch( void ) const
    {

         // Declare structure to be returned
        gnssDataMap toReturn;

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get the 'CommonTime' of the first element
            CommonTime firstEpoch( (*(*this).begin()).first );

            // Find the position of the first element PAST
            // 'firstEpoch+tolerance'
            gnssDataMap::const_iterator endPos(
                                 (*this).upper_bound(firstEpoch+tolerance) );

            // Add values to 'toReturn'
            for(gnssDataMap::const_iterator it = (*this).begin();
                it != endPos;
                ++it)
            {
                toReturn.insert( (*it) );
            }

        }

        return toReturn;

    }  // End of method 'gnssDataMap::frontEpoch()'



      // Removes the first epoch in the map. Be aware that this method
      // takes into account 'tolerance', so more than just one epoch may be
      // removed if they are within 'tolerance' margin.
    void gnssDataMap::pop_front_epoch( void )
    {

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get the 'CommonTime' of the first element
            CommonTime firstEpoch( (*(*this).begin()).first );

            // Find the position of the first element PAST
            // 'firstEpoch+tolerance'
            gnssDataMap::iterator endPos(
                                 (*this).upper_bound(firstEpoch+tolerance) );

            // Remove values
            for(gnssDataMap::iterator pos = (*this).begin();
                pos != endPos; )
            {

                // It is advisable to avoid sawing off the branch we are
                // sitting on
                (*this).erase( pos++ );
            }

        }

        return;

    }  // End of method 'gnssDataMap::pop_front_epoch()'



      // Return a copy of the last element in the map.
    gnssDataMap gnssDataMap::back( void ) const
    {

         // Object to be returned
        gnssDataMap toReturn;

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get last element
            gnssDataMap::const_reverse_iterator lastIt( (*this).rbegin() );

            // Insert the last element in 'toReturn' object
            toReturn.insert(*lastIt);

        }

        return toReturn;

    }  // End of method 'gnssDataMap::back()'



      // Removes the last element in the map.
    void gnssDataMap::pop_back( void )
    {

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Define an iterator pointing to the end of this map
            gnssDataMap::iterator pos( (*this).end() );

            // Delete element
            (*this).erase( --pos );

        }

        return;

    }  // End of method 'gnssDataMap::pop_back()'



      // Return the data corresponding to the last epoch in the map,
      // taking into account the 'tolerance' value.
    gnssDataMap gnssDataMap::backEpoch( void ) const
    {

         // Declare structure to be returned
        gnssDataMap toReturn;

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get the 'CommonTime' of the last element
            CommonTime lastEpoch( (*(--(*this).end())).first );

            // Find the position of the last element PAST
            // 'lastEpoch-tolerance'
            gnssDataMap::const_iterator beginPos(
                                 (*this).lower_bound(lastEpoch-tolerance) );

            // Add values to 'toReturn'
            for(gnssDataMap::const_iterator it = beginPos;
                it != (*this).end();
                ++it)
            {
                toReturn.insert( (*it) );
            }

        }

        return toReturn;

    }  // End of method 'gnssDataMap::backEpoch()'



      // Removes the last epoch in the map. Be aware that this method
      // takes into account 'tolerance', so more than just one epoch may be
      // removed if they are within 'tolerance' margin.
    void gnssDataMap::pop_back_epoch( void )
    {

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Get the 'CommonTime' of the last element
            CommonTime lastEpoch( (*(--(*this).end())).first );

            // Find the position of the last element PAST
            // 'lastEpoch-tolerance'
            gnssDataMap::iterator beginPos(
                                 (*this).lower_bound(lastEpoch-tolerance) );

            // Remove values
            for(gnssDataMap::iterator pos = beginPos;
                pos != (*this).end(); )
            {

                // It is advisable to avoid sawing off the branch we are
                // sitting on
                (*this).erase( pos++ );
            }

        }

        return;

    }  // End of method 'gnssDataMap::pop_back_epoch()'



      /* Return a 'gnssDataMap' with the data corresponding to provided
       * CommonTime, taking into account 'tolerance'.
       *
       * @param epoch         Epoch to be looked for.
       */
    gnssDataMap gnssDataMap::getDataFromEpoch( const CommonTime& epoch ) const
        throw( CommonTimeNotFound )
    {

         // Declare structure to be returned
        gnssDataMap toReturn;

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Find the position of the first element whose key is not less than
            // 'epoch-tolerance'
            gnssDataMap::const_iterator beginPos(
                                 (*this).lower_bound(epoch-tolerance) );

            // Find the position of the first element PAST
            // 'epoch+tolerance'
            gnssDataMap::const_iterator endPos(
                                 (*this).upper_bound(epoch+tolerance) );

            // Add values to 'toReturn'
            for(gnssDataMap::const_iterator it = beginPos;
                it != endPos;
                ++it)
            {
                toReturn.insert( (*it) );
            }

        }
        else
        {
            GPSTK_THROW(CommonTimeNotFound("Data map is empty"));
        }

         // Check if 'toReturn' is empty
        if( toReturn.empty() )
        {
            GPSTK_THROW(CommonTimeNotFound("Epoch not found"));
        }

        return toReturn;

    }  // End of method 'gnssDataMap::getDataFromEpoch()'



      /* Return the data value (double) corresponding to provided CommonTime,
       * SourceID, SatID and TypeID.
       *
       * @param epoch         Epoch to be looked for.
       * @param source        Source to be looked for.
       * @param satellite     Satellite to be looked for.
       * @param type          Type to be looked for.
       *
       * @warning If within (epoch +/- tolerance) more than one match exists,
       * then only the first one is returned.
       */
    double gnssDataMap::getValue( const CommonTime& epoch,
                                  const SourceID& source,
                                  const SatID& satellite,
                                  const TypeID& type ) const
        throw( CommonTimeNotFound, ValueNotFound )
    {

         // Look for the epoch (CommonTime) data
        gnssDataMap gdMap( getDataFromEpoch(epoch) );

         // Value to be returned
        double toReturn;

         // We'll need a flag
        bool found(false);

         // Look into the data structure
        for(gnssDataMap::const_iterator it = gdMap.begin();
            it != gdMap.end() && !found;
            ++it)
        {

            try
            {
                toReturn = (*it).second.getValue( source, satellite, type );
                found = true;
            }
            catch(...)
            {
                continue;
            }

        }

         // Check if value was found
        if( !found )
        {
            GPSTK_THROW(ValueNotFound("Value not found"));
        }

        return toReturn;

    }  // End of method 'gnssDataMap::getValue()'



      /* Return the data value (double) corresponding to the first epoch
       * in the data structure, given SourceID, SatID and TypeID.
       *
       * @param source        Source to be looked for.
       * @param satellite     Satellite to be looked for.
       * @param type          Type to be looked for.
       *
       * @warning If within first epoch (epoch +/- tolerance) more than one
       * match exists, then only the first one is returned.
       */
    double gnssDataMap::getValue( const SourceID& source,
                                  const SatID& satellite,
                                  const TypeID& type ) const
        throw( ValueNotFound )
    {

         // Look for the epoch (CommonTime) data
        gnssDataMap gdMap( frontEpoch() );

         // Value to be returned
        double toReturn;

         // We'll need a flag
        bool found(false);

         // Look into the data structure
        for(gnssDataMap::const_iterator it = gdMap.begin();
            it != gdMap.end() && !found;
            ++it)
        {

            try
            {
                toReturn = (*it).second.getValue( source, satellite, type );
                found = true;
            }
            catch(...)
            {
                continue;
            }

        }

         // Check if value was found
        if( !found )
        {
            GPSTK_THROW(ValueNotFound("Value not found"));
        }

        return toReturn;

    }  // End of method 'gnssDataMap::getValue()'



      /* Inserts a data value (double) at the provided CommonTime, SourceID,
       * SatID and TypeID, taking into account 'tolerance'.
       *
       * @param epoch         Epoch to be looked for.
       * @param source        Source to be looked for.
       * @param satellite     Satellite to be looked for.
       * @param type          Type of the new data.
       * @param value         Value to be inserted.
       */
    gnssDataMap& gnssDataMap::insertValue( const CommonTime& epoch,
                                           const SourceID& source,
                                           const SatID& satellite,
                                           const TypeID& type,
                                           double value )
        throw( CommonTimeNotFound, ValueNotFound )
    {

         // First check that structure isn't empty
        if( !( (*this).empty() ) )
        {

            // Find the position of the first element whose key is not less than
            // 'epoch-tolerance'
            gnssDataMap::iterator beginPos( (*this).lower_bound(epoch-tolerance) );

            // Find the position of the first element PAST
            // 'epoch+tolerance'
            gnssDataMap::iterator endPos( (*this).upper_bound(epoch+tolerance) );

            // Check if we found some match within current tolerance
            if( beginPos != endPos )
            {

                // We'll need a flag
                bool found(false);

                // Look into the range
                for(gnssDataMap::iterator it = beginPos;
                    it != endPos && !found;
                    ++it)
                {

                    // Look for the source
                    sourceDataMap::iterator it2( (*it).second.find(source) );

                    if( it2 != (*it).second.end() )
                    {

                        // Look for the satellite
                        satTypeValueMap::iterator it3((*it2).second.find(satellite));

                        if( it3 != (*it2).second.end() )
                        {

                            // Insert type and value
                            (*it3).second[ type ] = value;

                            // Work is done, let's get out
                            found = true;

                        }  // End of 'if( it3 != (*it2).second.end() )'

                    }  // End of 'if( it2 != (*it).second.end() )'

                }  // End of 'for( gnssDataMap::iterator it = beginPos ...'

                // Check if we found a proper place to insert value
                if( !found )
                {
                    GPSTK_THROW( ValueNotFound("No proper place to insert value"));
                }

            }
            else
            {
                // No match found for CommonTime with current tolerance
                GPSTK_THROW(CommonTimeNotFound("Epoch not found within tolerance"));
            }

        }
        else
        {
            GPSTK_THROW(CommonTimeNotFound("Data map is empty"));
        }

        return (*this);

    }  // End of method 'gnssDataMap::insertValue()'



      /* Inserts a data value (double) in the first epoch of the data
       * structure with the given SourceID, SatID and TypeID.
       *
       * @param source        Source to be looked for.
       * @param satellite     Satellite to be looked for.
       * @param type          Type of the new data.
       * @param value         Value to be inserted.
       */
    gnssDataMap& gnssDataMap::insertValue( const SourceID& source,
                                           const SatID& satellite,
                                           const TypeID& type,
                                           double value )
        throw( ValueNotFound )
    {

         // We'll need a flag
        bool found(false);

         // Look into the data structure
        for(gnssDataMap::iterator it = (*this).begin();
            it != (*this).end() && !found;
            ++it)
        {
            // Look for the source
            sourceDataMap::iterator it2( (*it).second.find(source) );

            if( it2 != (*it).second.end() )
            {
                // Look for the satellite
                satTypeValueMap::iterator it3( (*it2).second.find(satellite) );

                if( it3 != (*it2).second.end() )
                {
                    // Insert type and value
                    (*it3).second[ type ] = value;

                    // Work is done, let's get out
                    found = true;

                }  // End of 'if( it3 != (*it2).second.end() )'

            }  // End of 'if( it2 != (*it).second.end() )'

        }  // End of 'for( gnssDataMap::iterator it = (*this).begin(); ...'

         // Check if we found a proper place to insert value
        if( !found )
        {
            GPSTK_THROW( ValueNotFound("No proper place to insert value"));
        }

        return (*this);

    }  // End of method 'gnssDataMap::insertValue()'



      /* Get a set with all the SourceID's in this data structure.
       *
       * @warning If current 'gnssDataMap' is big, this could be a very
       * costly operation.
       */
    SourceIDSet gnssDataMap::getSourceIDSet( void ) const
    {

         // SourceID set to be returned
        SourceIDSet toReturn;

         // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = (*this).begin();
            it != (*this).end();
            ++it)
        {

            // Then, iterate through corresponding 'sourceDataMap'
            for(sourceDataMap::const_iterator sdmIter = (*it).second.begin();
                sdmIter != (*it).second.end();
                ++sdmIter)
            {

                // Add current SourceID to 'toReturn'
                toReturn.insert( (*sdmIter).first );

            }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return toReturn;

    }  // End of method 'gnssDataMap::getSourceIDSet()'


    SourceIDSet gnssDataMap::getSourceIDSet( const SatID& satellite,
                                             const CommonTime& time )
        throw( CommonTimeNotFound )
    {

        gnssDataMap gdsMap( getDataFromEpoch( time ) );

        SourceIDSet  sourceIDSet;

        if( !gdsMap.empty() )
        {

            for(gnssDataMap::const_iterator gdsMapIt = gdsMap.begin();
                gdsMapIt != gdsMap.end();
                ++gdsMapIt)
            {

                for(sourceDataMap::const_iterator sdMapIt = gdsMapIt->second.begin();
                    sdMapIt != gdsMapIt->second.end();
                    ++sdMapIt)
                {

                    /// find satellite in map
                    satTypeValueMap::const_iterator stvMapIt = sdMapIt->second.find( satellite );

                    if( stvMapIt != sdMapIt->second.end() )
                    {

                        /// insert SourceID to Set
                        sourceIDSet.insert( sdMapIt->first );
                        break;

                    }

                } // End of "sourceDataMap::const_iterator sdMapIt = ..."

            } // End of "gnssDataMap::const_iterator gdsMapIt = ..."

        }

        return sourceIDSet;

    }


      /* Get a set with all the SatID's in this data structure.
       *
       * @warning If current 'gnssDataMap' is big, this could be a very
       * costly operation.
       */
    SatIDSet gnssDataMap::getSatIDSet( void ) const
    {

         // SatID set to be returned
        SatIDSet toReturn;

         // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = (*this).begin();
            it != (*this).end();
            ++it)
        {

            // Then, iterate through corresponding 'sourceDataMap'
            for(sourceDataMap::const_iterator sdmIter = (*it).second.begin();
                sdmIter != (*it).second.end();
                ++sdmIter)
            {

                // Finally, iterate through corresponding 'satTypeValueMap'
                for(satTypeValueMap::const_iterator stvmIter =
                                                      (*sdmIter).second.begin();
                    stvmIter != (*sdmIter).second.end();
                    ++stvmIter)
                {

                    // Add current SatID to 'toReturn'
                    toReturn.insert( (*stvmIter).first );

                }  // End of 'for( satTypeValueMap::const_iterator stvmIter = ...'

            }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return toReturn;

    }  // End of method 'gnssDataMap::getSatIDSet()'



      // Convenience output method
    std::ostream& gnssDataMap::dump( std::ostream& s, int mode ) const
    {

         // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = (*this).begin();
            it != (*this).end();
            ++it)
        {

            // Then, iterate through corresponding 'sourceDataMap'
            for(sourceDataMap::const_iterator sdmIter = (*it).second.begin();
                sdmIter != (*it).second.end();
                ++sdmIter)
            {

                // Finally, iterate through corresponding 'satTypeValueMap'
                for(satTypeValueMap::const_iterator stvmIter =
                                                      (*sdmIter).second.begin();
                    stvmIter != (*sdmIter).second.end();
                    ++stvmIter)
                {

                    // Declare a 'YDSTime' object to ease printing
                    YDSTime time( (*it).first );

                    // First, print year, Day-Of-Year and Seconds of Day
                    s << time.year << " "
                      << time.doy << " "
                      << time.sod << " ";

                    // Second, print SourceID information
                    s << (*sdmIter).first << " ";

                    // Third, print satellite (system and PRN)
                    s << (*stvmIter).first << " ";

                    // Fourth, print type descriptions and numerical values
                    for(typeValueMap::const_iterator itObs =
                                                   (*stvmIter).second.begin();
                        itObs != (*stvmIter).second.end();
                        ++itObs)
                    {

                        s << (*itObs).first << " "
                          << (*itObs).second << " ";

                    }  // End of 'for( typeValueMap::const_iterator itObs = ...'

                    // Add end-of-line
                    s << endl;

                }  // End of 'for( satTypeValueMap::const_iterator stvmIter = ...'

            }  // End of 'for( sourceDataMap::const_iterator sdmIter = ...'

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

         // Let's return the 'std::ostream'
        return s;

    }  // End of method 'gnssDataMap::dump()'



      // Return a gnssDataMap with only this source.
      // @param source Source to be extracted.
    gnssDataMap gnssDataMap::extractSourceID(const SourceID& source)
    {
        SourceIDSet sourceSet;
        sourceSet.insert(source);

        return extractSourceID(sourceSet);

    }  // End of method 'gnssDataMap::extractSourceID()'



      /// Return a gnssDataMap with only these sources.
      /// @param sourceSet Set(SourceIDSet) containing the sources
      ///                  to be extracted.
    gnssDataMap gnssDataMap::extractSourceID(const SourceIDSet& sourceSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                SourceIDSet::const_iterator itsrc2 = sourceSet.find(itsrc->first);
                if(itsrc2!=sourceSet.end())
                {
                    gnssSatTypeValue gds;
                    gds.header.epoch = time;
                    gds.header.source = itsrc->first;
                    gds.body = itsrc->second;

                    dataMap.addGnssSatTypeValue(gds);
                }

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return dataMap;

    }  // End of method 'gnssDataMap::extractSourceID()'


      // Modifies this object, keeping only this source.
      // @param source Source to be extracted.
    gnssDataMap& gnssDataMap::keepOnlySourceID(const SourceID& source)
    {
        (*this) = extractSourceID(source);
        return (*this);
    }


      // Modifies this object, keeping only these sources.
      // @param sourceSet Set(SourceIDSet) containing the sources
      //                  to be extracted.
    gnssDataMap& gnssDataMap::keepOnlySourceID(const SourceIDSet& sourceSet)
    {
        (*this) = extractSourceID(sourceSet);
        return (*this);
    }


      // Modifies this object, removing this source.
      // @param source Source to be removed.
    gnssDataMap& gnssDataMap::removeSourceID(const SourceID& source)
    {
        SourceIDSet sourceSet;
        sourceSet.insert(source);

        return removeSourceID(sourceSet);
    }


      // Modifies this object, keeping only these sources.
      // @param sourceSet Set(SourceIDSet) containing the sources
      //                  to be removed.
    gnssDataMap& gnssDataMap::removeSourceID(const SourceIDSet& sourceSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                SourceIDSet::const_iterator itsrc2 = sourceSet.find(itsrc->first);
                if(itsrc2==sourceSet.end())
                {
                    gnssSatTypeValue gds;
                    gds.header.epoch = time;
                    gds.header.source = itsrc->first;
                    gds.body = itsrc->second;

                    dataMap.addGnssSatTypeValue(gds);
                }

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        (*this) = dataMap;

        return (*this);

    }  // End of method 'gnssDataMap::removeSourceID()'


      // Return a gnssDataMap with only this satellite.
      // @param sat Satellite to be extracted.
    gnssDataMap gnssDataMap::extractSatID(const SatID& sat)
    {
        SatIDSet satSet;
        satSet.insert(sat);

        return extractSatID(satSet);
    }


      // Return a gnssDataMap with only these satellites.
      // @param satSet Set(SatIDSet) containing the satellite
      //               to be extracted.
    gnssDataMap gnssDataMap::extractSatID(const SatIDSet& satSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue gds;
                gds.header.epoch = time;
                gds.header.source = itsrc->first;
                gds.body = itsrc->second;

                gds.body.keepOnlySatID(satSet);

                dataMap.addGnssSatTypeValue(gds);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return dataMap;

    }  // End of method 'gnssDataMap::extractSatID()'


      // Modifies this object, keeping only this satellite.
      // @param sat Satellite to be extracted.
    gnssDataMap& gnssDataMap::keepOnlySatID(const SatID& sat)
    {
        (*this) = extractSatID(sat);
        return (*this);
    }


      // Modifies this object, keeping only these satellites.
      // @param satSet Set(SatIDSet) containing the satellite
      //                  to be extracted.
    gnssDataMap& gnssDataMap::keepOnlySatID(const SatIDSet& satSet)
    {
        (*this) = extractSatID(satSet);
        return (*this);
    }


      // Modifies this object, removing this satellite.
      // @param sat Satellite to be removed.
    gnssDataMap& gnssDataMap::removeSatID(const SatID& sat)
    {
        SatIDSet satSet;
        satSet.insert(sat);

        return removeSatID(satSet);
    }


      // Modifies this object, keeping only these satellites.
      // @param satSet Set(SatIDSet) containing the satellites
      //               to be removed.
    gnssDataMap& gnssDataMap::removeSatID(const SatIDSet& satSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue gds;
                gds.header.epoch = time;
                gds.header.source = itsrc->first;
                gds.body = itsrc->second;

                gds.body.removeSatID(satSet);

                dataMap.addGnssSatTypeValue(gds);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        (*this) = dataMap;

        return (*this);

    }  // End of method 'gnssDataMap::removeSatID()'


      // Modifies this object, removing this satellite from this source.
      // @param source Source to be removed.
      // @param sat Satellite to be removed.
    gnssDataMap& gnssDataMap::removeSatID( const SourceID& source,
                                           const SatID& sat )
    {
        for(gnssDataMap::iterator iter = this->begin();
            iter != this->end();
            ++iter)
        {
            sourceDataMap::iterator sdmIter = iter->second.find( source );

            if( iter->second.end() != sdmIter )
            {
                sdmIter->second.removeSatID(sat);
            }
        }

    }  // End of method 'gnssDataMap::removeSatID()'


      /// Return a gnssDataMap with only this type.
      /// @param type Type to be extracted.
    gnssDataMap gnssDataMap::extractTypeID(const TypeID& type)
    {
        TypeIDSet typeSet;
        typeSet.insert(type);

        return extractTypeID(typeSet);
    }


      // Return a gnssDataMap with only these satellites.
      // @param typeSet Set(TypeIDSet) containing the types
      //               to be extracted.
    gnssDataMap gnssDataMap::extractTypeID(const TypeIDSet& typeSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue gds;
                gds.header.epoch = time;
                gds.header.source = itsrc->first;
                gds.body = itsrc->second;

                gds.body.keepOnlyTypeID(typeSet);

                dataMap.addGnssSatTypeValue(gds);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        return dataMap;
    }


      // Modifies this object, keeping only this type.
      // @param type Type to be extracted.
    gnssDataMap& gnssDataMap::keepOnlyTypeID(const TypeID& type)
    {
        (*this) = extractTypeID(type);
        return (*this);
    }


      // Modifies this object, keeping only these types.
      // @param typeSet Set(TypeIDSet) containing the type
      //                  to be extracted.
    gnssDataMap& gnssDataMap::keepOnlyTypeID(const TypeIDSet& typeSet)
    {
        (*this) = extractTypeID(typeSet);
        return (*this);
    }


      // Modifies this object, removing this type.
      // @param type Type to be removed.
    gnssDataMap& gnssDataMap::removeTypeID(const TypeID& type)
    {
        TypeIDSet typeSet;
        typeSet.insert(type);

        return removeTypeID(typeSet);
    }


      // Modifies this object, remove these types.
      // @param typeSet Set(TypeIDSet) containing the types
      //               to be removed.
    gnssDataMap& gnssDataMap::removeTypeID(const TypeIDSet& typeSet)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue gds;
                gds.header.epoch = time;
                gds.header.source = itsrc->first;
                gds.body = itsrc->second;

                gds.body.removeTypeID(typeSet);

                dataMap.addGnssSatTypeValue(gds);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        (*this) = dataMap;

        return (*this);

    }  // End of method 'gnssDataMap::removeTypeID()'


     /// Modifies this object, remove this satellite of this source
     /// and this satellite.
     /// @param source  Source to be removed.
     /// @param sat     Satellite to be removed.
     /// @param type    Type to be removed.
    gnssDataMap& gnssDataMap::removeTypeID(const SourceID& source,
                                           const SatID& sat,
                                           const TypeID& type)
    {
        gnssDataMap dataMap;

        // Iterate through all items in the gnssDataMap
        for(gnssDataMap::const_iterator it = this->begin();
            it != this->end();
            ++it)
        {
            const CommonTime& time(it->first);
            const sourceDataMap& sourceMap(it->second);

            for(sourceDataMap::const_iterator itsrc = sourceMap.begin();
                itsrc != sourceMap.end();
                ++itsrc)
            {
                gnssSatTypeValue gds;
                gds.header.epoch = time;
                gds.header.source = itsrc->first;
                gds.body = itsrc->second;

                if(itsrc->first == source)
                {
                    for(satTypeValueMap::iterator itsat = gds.body.begin();
                        itsat != gds.body.end();
                        ++itsat)
                    {
                        if(itsat->first == sat)
                        {
                            itsat->second.erase(type);
                        }
                    }
                }

                dataMap.addGnssSatTypeValue(gds);

            }  // loop in the sources

        }  // End of 'for( gnssDataMap::const_iterator it = gdMap.begin(); ...'

        (*this) = dataMap;

        return (*this);

    }  // End of method 'gnssDataMap::removeTypeID()'


      /* Edit the dataset, removing data outside the indicated time
       *  interval.
       *
       * @param[in] tmin defines the beginning of the time interval
       * @param[in] tmax defines the end of the time interval
       */
    gnssDataMap& gnssDataMap::edit(CommonTime tmin, CommonTime tmax )
    {
        gnssDataMap dataMap;

        while(!this->empty())
        {
            gnssDataMap gds = this->frontEpoch();

            CommonTime time(gds.begin()->first);
            if( (time>=tmin) && (time<=tmax) ) dataMap.addGnssDataMap(gds);

            this->pop_front_epoch();
        }

        (*this) = dataMap;

        return (*this);
    }



      // Load data from a rinex observation file
    void gnssDataMap::loadObsFile(std::string obsFile)
    {
        try
        {
            RinexObsStream rin(obsFile);
            rin.exceptions(ios::failbit);

            gnssRinex gRin;
            while( rin >> gRin )
            {
                this->addGnssRinex(gRin);
            }

            rin.close();
        }
        catch(...)
        {
            Exception e("Failed to load obs file '"+obsFile+"'.");
            GPSTK_THROW(e);
        }

    }  // End of method 'gnssDataMap::loadObsFile()'



      ////// Other stuff //////



      // Convenience output method for structure satTypeValueMap
    std::ostream& satTypeValueMap::dump( std::ostream& s, int mode ) const
    {

        for(satTypeValueMap::const_iterator it = (*this).begin();
            it != (*this).end();
            ++it)
        {

            // First, print satellite (system and PRN)
            s << (*it).first << " ";

            for(typeValueMap::const_iterator itObs = (*it).second.begin();
                itObs != (*it).second.end();
                ++itObs)
            {

                if (mode==1)
                {
                    s << (*itObs).first << " ";
                }

                s << (*itObs).second << " ";

            }  // End of 'for( typeValueMap::const_iterator itObs = ...'

            s << endl;

        }  // End of 'for( satTypeValueMap::const_iterator it = ...'

         // Let's return the 'std::ostream'
        return s;

    }  // End of method 'satTypeValueMap::dump()'



      // stream output for satTypeValueMap
    std::ostream& operator<<( std::ostream& s, const satTypeValueMap& stvMap )
    {

        stvMap.dump(s);
        return s;

    }  // End of 'operator<<'



      // stream output for gnssDataMap
    std::ostream& operator<<( std::ostream& s, const gnssDataMap& gdsMap )
    {

        gdsMap.dump(s);
        return s;

    }  // End of 'operator<<'


      // Stream input for gnssRinex
    std::istream& operator>>( std::istream& i, gnssRinex& f )
    {
        if( RinexObsStream::isRinexObsStream(i) )       // Rinex2
        {
            try
            {
                RinexObsStream& strm = dynamic_cast<RinexObsStream&>(i);

                if(!strm.headerRead) strm >> strm.header;

                RinexObsHeader& roh = strm.header;

                RinexObsData rod;
                strm >> rod;

                f.header.source.type = SatIDsystem2SourceIDtype(roh.system);
                f.header.source.sourceName = roh.markerName.substr(0,4);

                f.header.antennaType = roh.antType;
                f.header.antennaPosition = roh.antennaPosition;
                f.header.epochFlag = rod.epochFlag;
                f.header.epoch = rod.time;

                f.body = satTypeValueMapFromRinexObsData(roh, rod);

                return i;
            }
            catch(...)
            {
                return i;
            }
        }

        if( Rinex3ObsStream::isRinex3ObsStream(i) )     // Rinex3
        {
            Rinex3ObsStream& strm = dynamic_cast<Rinex3ObsStream&>(i);

            // If the header hasn't been read, read it...
            if(!strm.headerRead) strm >> strm.header;

            // Clear out this object
            Rinex3ObsHeader& roh = strm.header;

            Rinex3ObsData rod;
            strm >> rod;

            // Fill data
            f.header.source.type = SatIDsystem2SourceIDtype(roh.fileSysSat);
            f.header.source.sourceName = roh.markerName.substr(0,4);
//            f.header.source.sourceNumber = roh.markerNumber;

            f.header.antennaType = roh.antType;
            f.header.antennaPosition = roh.antennaPosition;
            f.header.epochFlag = rod.epochFlag;
            f.header.epoch = rod.time;

            f.body = satTypeValueMapFromRinex3ObsData(roh, rod);

            return i;
        }

    }  // End of stream input for gnssRinex



      // Convenience function to convert from SatID system to SourceID type.
      // @param sid Satellite ID.
    SourceID::SourceType SatIDsystem2SourceIDtype(const SatID& sid)
    {

         // Select the right system the data came from
        switch(sid.system)
        {
            case SatID::systemGPS:
                return SourceID::GPS;
                break;
            case SatID::systemGLONASS:
                return SourceID::GLONASS;
                break;
            case SatID::systemGalileo:
                return SourceID::Galileo;
                break;
            case SatID::systemSBAS:
                return SourceID::SBAS;
                break;
            case SatID::systemQZSS:
                return SourceID::QZSS;
                break;
            case SatID::systemBDS:
                return SourceID::BDS;
                break;
            case SatID::systemIRNSS:
                return SourceID::IRNSS;
                break;
            case SatID::systemMixed:
                return SourceID::Mixed;
                break;
            default:
                return SourceID::Unknown;
        }

    } // End SatIDsystem2SourceIDtype(const SatID& sid)


      // Convenience function to fill a satTypeValueMap with data
      // from RinexObsData.
      // @param roh RinexObsHeader holding the data
      // @param rod RinexObsData holding the data.
    satTypeValueMap satTypeValueMapFromRinexObsData(const RinexObsHeader& roh,
                                                    const RinexObsData& rod)
    {
        // We need to declare a satTypeValueMap
        satTypeValueMap theMap;

        RinexObsData::RinexSatMap::const_iterator it;
        for(RinexObsData::RinexSatMap::const_iterator it = rod.obs.begin();
            it != rod.obs.end();
            ++it)
        {
            RinexObsData::RinexObsTypeMap otmap = (*it).second;
            SatID sat = (*it).first;

            typeValueMap tvMap;

            for(RinexObsData::RinexObsTypeMap::const_iterator itObs = otmap.begin();
                itObs != otmap.end();
                ++itObs)
            {
                RinexSatID rsat(sat.id,sat.system);

                TypeID type = ConvertToTypeID( itObs->first, rsat );

                const bool isPhase = IsCarrierPhase( itObs->first );
                const int n = GetCarrierBand( itObs->first );

                if(isPhase)
                {
                    if( (*itObs).second.data == 0.0 ) continue;

                    tvMap[ type ] = (*itObs).second.data*getWavelength(rsat,n);
//                    tvMap[ type ] = (*itObs).second.data;     // for debug

                    if(n == 1)
                    {
                        tvMap[TypeID::LLI1] = (*itObs).second.lli;
                        tvMap[TypeID::SSI1] = (*itObs).second.ssi;
                    }
                    else if(n == 2)
                    {
                        tvMap[TypeID::LLI2] = (*itObs).second.lli;
                        tvMap[TypeID::SSI2] = (*itObs).second.ssi;
                    }
                    else if(n == 3)
                    {
                        tvMap[TypeID::LLI3] = (*itObs).second.lli;
                        tvMap[TypeID::SSI3] = (*itObs).second.ssi;
                    }
                    else if(n == 5)
                    {
                        tvMap[TypeID::LLI5] = (*itObs).second.lli;
                        tvMap[TypeID::SSI5] = (*itObs).second.ssi;
                    }
                    else if(n == 6)
                    {
                        tvMap[TypeID::LLI6] = (*itObs).second.lli;
                        tvMap[TypeID::SSI6] = (*itObs).second.ssi;
                    }
                    else if(n == 7)
                    {
                        tvMap[TypeID::LLI7] = (*itObs).second.lli;
                        tvMap[TypeID::SSI7] = (*itObs).second.ssi;
                    }
                    else if(n == 8)
                    {
                        tvMap[TypeID::LLI8] = (*itObs).second.lli;
                        tvMap[TypeID::SSI8] = (*itObs).second.ssi;
                    }
                }
                else
                {
                    if( (*itObs).second.data == 0.0 ) continue;

                    tvMap[ type ] = (*itObs).second.data;
                }
            }

            theMap[sat] = tvMap;

        }   // End loop over all the satellite

        return theMap;

    }


      // Convenience function to fill a satTypeValueMap with data
      // from Rinex3ObsData.
      // @param roh Rinex3ObsHeader holding the data
      // @param rod Rinex3ObsData holding the data.
    satTypeValueMap satTypeValueMapFromRinex3ObsData(const Rinex3ObsHeader& roh,
                                                     const Rinex3ObsData& rod)
    {

        // We need to declare a satTypeValueMap
        satTypeValueMap theMap;

        for(Rinex3ObsData::DataMap::const_iterator it = rod.obs.begin();
            it != rod.obs.end();
            ++it)
        {
            RinexSatID sat(it->first);

            typeValueMap tvMap;

            map<std::string,std::vector<RinexObsID> > mapObsTypes(roh.mapObsTypes);
            const vector<RinexObsID> types = mapObsTypes[sat.toString().substr(0,1)];

            for(size_t i=0; i<types.size(); i++)
            {
                TypeID type = ConvertToTypeID(types[i],sat);

                const int n = GetCarrierBand(types[i]);

                if(types[i].type == ObsID::otC)         // Range
                {
                    if(it->second[i].data == 0.0) continue;

                    tvMap[ type ] = it->second[i].data;
                }
                else if(types[i].type == ObsID::otL)    // Phase
                {
                    if(it->second[i].data == 0.0) continue;

                    tvMap[type] = it->second[i].data*getWavelength(sat,n);

                    // n=1 2 5 6 7 8 9
                    if(n == 1)
                    {
                        tvMap[TypeID::LLI1] = it->second[i].lli;
                        tvMap[TypeID::SSI1] = it->second[i].ssi;
                    }
                    else if(n == 2)
                    {
                        tvMap[TypeID::LLI2] = it->second[i].lli;
                        tvMap[TypeID::SSI2] = it->second[i].ssi;
                    }
                    else if(n == 3)
                    {
                        tvMap[TypeID::LLI3] = it->second[i].lli;
                        tvMap[TypeID::SSI3] = it->second[i].ssi;
                    }
                    else if(n == 5)
                    {
                        tvMap[TypeID::LLI5] = it->second[i].lli;
                        tvMap[TypeID::SSI5] = it->second[i].ssi;
                    }
                    else if(n == 6)
                    {
                        tvMap[TypeID::LLI6] = it->second[i].lli;
                        tvMap[TypeID::SSI6] = it->second[i].ssi;
                    }
                    else if(n == 7)
                    {
                        tvMap[TypeID::LLI7] = it->second[i].lli;
                        tvMap[TypeID::SSI7] = it->second[i].ssi;
                    }
                    else if(n == 8)
                    {
                        tvMap[TypeID::LLI8] = it->second[i].lli;
                        tvMap[TypeID::SSI8] = it->second[i].ssi;
                    }
                    else if(n == 9)
                    {
                        tvMap[TypeID::LLI9] = it->second[i].lli;
                        tvMap[TypeID::SSI9] = it->second[i].ssi;
                    }
                }
                else
                {
                    tvMap[ type ] = it->second[i].data;
                }
            }

            theMap[sat] = tvMap;

        }   // End loop over all the satellite

        return theMap;

    }


//    satTypeValueMap satTypeValueMapFromRinex3ObsData(
//                         const Rinex3ObsHeader& roh, const Rinex3ObsData& rod )
//    {
//        // Valid Rinex Tracking Codes of GPS
//        const string G1_validRTCs( ObsID::validRinexTrackingCodes['G']['1'] );
//        const string G2_validRTCs( ObsID::validRinexTrackingCodes['G']['2'] );
//        const string G5_validRTCs( ObsID::validRinexTrackingCodes['G']['5'] );
//
//        // Valid Rinex Tracking Codes of Galileo
//        const string E1_validRTCs( ObsID::validRinexTrackingCodes['E']['1'] );
//        const string E5_validRTCs( ObsID::validRinexTrackingCodes['E']['5'] );
//        const string E6_validRTCs( ObsID::validRinexTrackingCodes['E']['6'] );
//        const string E7_validRTCs( ObsID::validRinexTrackingCodes['E']['7'] );
//        const string E8_validRTCs( ObsID::validRinexTrackingCodes['E']['8'] );
//
//        // Valid Rinex Tracking Codes of BDS
//        const string C1_validRTCs( ObsID::validRinexTrackingCodes['C']['1'] );
//        const string C2_validRTCs( ObsID::validRinexTrackingCodes['C']['2'] );
//        const string C6_validRTCs( ObsID::validRinexTrackingCodes['C']['6'] );
//        const string C7_validRTCs( ObsID::validRinexTrackingCodes['C']['7'] );
//
//        // Tracking Codes in Observation File Header
//        vector<char> G1_obsTCs, G2_obsTCs, G5_obsTCs;
//        vector<char> E1_obsTCs, E5_obsTCs, E6_obsTCs, E7_obsTCs, E8_obsTCs;
//        vector<char> C2_obsTCs, C6_obsTCs, C7_obsTCs;
//
//        // Observation Types to be used
//        vector<RinexObsID> useObsTypesOfGPS;
//        vector<RinexObsID> useObsTypesOfGalileo;
//        vector<RinexObsID> useObsTypesOfBDS;
//
//        map< string,vector<RinexObsID> > mapObsTypes(roh.mapObsTypes);
//        map< string,vector<RinexObsID> >::const_iterator it_mot;
//
//        // Loop of systems
//        for(it_mot=mapObsTypes.begin(); it_mot!=mapObsTypes.end(); ++it_mot)
//        {
//            // system
//            string sys = it_mot->first;
//
//            // RinexObsID
//            vector<RinexObsID> roi = it_mot->second;
//
//            vector<RinexObsID>::const_iterator it_roi;
//
//            if(sys == "G") // GPS
//            {
//                // Loop of RinexObsID
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    // Count Trancking Codes Observed
//                    if(it_roi->type == ObsID::otC)
//                    {
//                        char tc;
//
//                        if(it_roi->band==ObsID::cb1)         // G1
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            G1_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb2)    // G2
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            G2_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb5)    // G5
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            G5_obsTCs.push_back( tc );
//                        }
//                    }
//                }
//
//                // Determine Tracking Code to be Used
//                char G1_useTC(' ');
//                if(G1_obsTCs.size() > 0)
//                {
//                    G1_useTC = G1_obsTCs[0];
//                    for(size_t i=0; i<G1_obsTCs.size(); i++)
//                    {
//                        if( G1_validRTCs.find(G1_obsTCs[i]) < G1_validRTCs.find(G1_useTC) )
//                            G1_useTC = G1_obsTCs[i];
//                    }
//                }
//
//                char G2_useTC(' ');
//                if(G2_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<G2_obsTCs.size(); i++)
//                    {
//                        if( G2_validRTCs.find(G2_obsTCs[i]) < G2_validRTCs.find(G2_useTC) )
//                           G2_useTC = G2_obsTCs[i];
//                    }
//                }
//
//                char G5_useTC(' ');
//                if(G5_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<G5_obsTCs.size(); i++)
//                    {
//                        if(G5_validRTCs.find(G5_obsTCs[i]) < G5_validRTCs.find(G5_useTC))
//                            G5_useTC = G5_obsTCs[i];
//                    }
//                }
//
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    if((it_roi->band==ObsID::cb1 &&
//                        it_roi->code==ObsID::char2tc.find(G1_useTC)->second)
//                    || (it_roi->band==ObsID::cb2 &&
//                        it_roi->code==ObsID::char2tc.find(G2_useTC)->second)
//                    || (it_roi->band==ObsID::cb5 &&
//                        it_roi->code==ObsID::char2tc.find(G5_useTC)->second))
//                    {
//                        useObsTypesOfGPS.push_back(*it_roi);
//                    }
//                }
////               cout << "G: ";
////               cout << G1_useTC << ' ' << G2_useTC << ' ' << G5_useTC << endl;
//            }
//            else if(sys == "E")
//            {
//                // Loop of RinexObsID
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    // Count Trancking Codes Observed
//                    if(it_roi->type == ObsID::otC)
//                    {
//                        char tc;
//
//                        if(it_roi->band==ObsID::cb1)         // E1
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            E1_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb5)    // E5
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            E5_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb6)    // E6
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            E6_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb7)    // E7
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            E7_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb8)    // E8
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            E8_obsTCs.push_back( tc );
//                        }
//                    }
//                }
//
//                // Determine Tracking Code to be Used
//
//                char E1_useTC(' ');
//                if(E1_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<E1_obsTCs.size(); i++)
//                    {
//                        if(E1_validRTCs.find(E1_obsTCs[i]) < E1_validRTCs.find(E1_useTC))
//                            E1_useTC = E1_obsTCs[i];
//                    }
//                }
//
//                char E5_useTC(' ');
//                if(E5_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<E5_obsTCs.size(); i++)
//                    {
//                        if(E5_validRTCs.find(E5_obsTCs[i]) < E5_validRTCs.find(E5_useTC))
//                            E5_useTC = E5_obsTCs[i];
//                    }
//                }
//
//                char E6_useTC(' ');
//                if(E6_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<E6_obsTCs.size(); i++)
//                    {
//                        if(E6_validRTCs.find(E6_obsTCs[i]) < E6_validRTCs.find(E6_useTC))
//                            E6_useTC = E6_obsTCs[i];
//                    }
//                }
//
//                char E7_useTC(' ');
//                if(E7_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<E7_obsTCs.size(); i++)
//                    {
//                        if(E7_validRTCs.find(E7_obsTCs[i]) < E7_validRTCs.find(E7_useTC))
//                            E7_useTC = E7_obsTCs[i];
//                    }
//                }
//
//                char E8_useTC(' ');
//                if(E8_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<E8_obsTCs.size(); i++)
//                    {
//                        if(E8_validRTCs.find(E8_obsTCs[i]) < E8_validRTCs.find(E8_useTC))
//                            E8_useTC = E8_obsTCs[i];
//                    }
//                }
//
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    if((it_roi->band == ObsID::cb1 &&
//                        it_roi->code == ObsID::char2tc.find(E1_useTC)->second)
//                    || (it_roi->band == ObsID::cb5 &&
//                        it_roi->code == ObsID::char2tc.find(E5_useTC)->second)
//                    || (it_roi->band == ObsID::cb6 &&
//                        it_roi->code == ObsID::char2tc.find(E6_useTC)->second)
//                    || (it_roi->band == ObsID::cb7 &&
//                        it_roi->code == ObsID::char2tc.find(E7_useTC)->second)
//                    || (it_roi->band == ObsID::cb8 &&
//                        it_roi->code == ObsID::char2tc.find(E8_useTC)->second))
//                    {
//                        useObsTypesOfGalileo.push_back(*it_roi);
//                    }
//                }
////               cout << "E: ";
////               cout << E1_useTC << ' ' << E5_useTC << ' ' << E6_useTC << ' '
////                    << E7_useTC << ' ' << E8_useTC << endl;
//            }
//            else if(sys == "C")
//            {
//                // Loop of RinexObsID
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    // Count Trancking Codes Observed
//                    if(it_roi->type == ObsID::otC)
//                    {
//                        char tc;
//
//                        if(it_roi->band==ObsID::cb1 ||
//                           it_roi->band==ObsID::cb2)         // B1
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            C2_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb6)    // B2
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            C6_obsTCs.push_back( tc );
//                        }
//                        else if(it_roi->band==ObsID::cb7)    // B3
//                        {
//                            tc = ObsID::tc2char.find(it_roi->code)->second;
//                            C7_obsTCs.push_back( tc );
//                        }
//                    }
//                }
//
//                // Determine Tracking Code to be Used
//
//                char C2_useTC(' ');
//                if(C2_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<C2_obsTCs.size(); i++)
//                    {
//                        if(C2_validRTCs.find(C2_obsTCs[i]) < C2_validRTCs.find(C2_useTC))
//                           C2_useTC = C2_obsTCs[i];
//                    }
//                }
//
//                char C6_useTC(' ');
//                if(C6_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<C6_obsTCs.size(); i++)
//                    {
//                        if(C6_validRTCs.find(C6_obsTCs[i]) < C6_validRTCs.find(C6_useTC))
//                            C6_useTC = C6_obsTCs[i];
//                    }
//                }
//
//                char C7_useTC(' ');
//                if(C7_obsTCs.size() > 0)
//                {
//                    for(size_t i=0; i<C7_obsTCs.size(); i++)
//                    {
//                        if(C7_validRTCs.find(C7_obsTCs[i]) < C7_validRTCs.find(C7_useTC))
//                            C7_useTC = C7_obsTCs[i];
//                    }
//                }
//
//                for(it_roi=roi.begin(); it_roi!=roi.end(); it_roi++)
//                {
//                    if(((it_roi->band==ObsID::cb1||it_roi->band==ObsID::cb2) &&
//                        it_roi->code==ObsID::char2tc.find(C2_useTC)->second)
//                    || (it_roi->band==ObsID::cb6 &&
//                        it_roi->code==ObsID::char2tc.find(C6_useTC)->second)
//                    || (it_roi->band==ObsID::cb7 &&
//                        it_roi->code==ObsID::char2tc.find(C7_useTC)->second))
//                    {
//                        useObsTypesOfBDS.push_back(*it_roi);
//                    }
//                }
////               cout << "C: ";
////               cout << C2_useTC << ' ' << C6_useTC << ' ' << C7_useTC << endl;
//            }
//
//        } // End of for(it_mot...)
//
//
//
//        // We need to declare a satTypeValueMap
//        satTypeValueMap theMap;
//
//        // map< RinexSatID, std::vector<RinexDatum> >
//        Rinex3ObsData::DataMap::const_iterator it;
//
//        for(it=rod.obs.begin(); it != rod.obs.end(); ++it)
//        {
//            RinexSatID sat(it->first);
//
//            typeValueMap tvMap;
//
//            vector<RinexObsID> roi_obs(mapObsTypes[sat.toString().substr(0,1)]);
//
//            vector<RinexObsID> roi_use;
//
//            if(sat.system == SatID::systemGPS)
//                roi_use = useObsTypesOfGPS;
//            else if(sat.system == SatID::systemGalileo)
//                roi_use = useObsTypesOfGalileo;
//            else if(sat.system == SatID::systemBDS)
//                roi_use = useObsTypesOfBDS;
//            else
//                continue;
//
////           cout << sat << ' ';
//
//            for(size_t i=0; i<roi_obs.size(); i++)
//            {
//                if( find(roi_use.begin(),roi_use.end(),roi_obs[i]) == roi_use.end() )
//                    continue;
//
//                TypeID type_3c( ConvertToTypeID(roi_obs[i],sat) );
//                TypeID type_2c( asString(type_3c).substr(0,2) );
//
////               cout << type_2c << ' ';
//
//                int n = GetCarrierBand(roi_obs[i]);
//                double data = it->second[i].data;
//
//                if(data == 0.0) continue;
//
//                if(roi_obs[i].type == ObsID::otL)
//                {
//                    tvMap[type_2c] = data * getWavelength(sat,n);
//
//                    double lli = it->second[i].lli;
//                    double ssi = it->second[i].ssi;
//
//                    if(n == 1)
//                    {
//                        tvMap[TypeID::LLI1] = lli; tvMap[TypeID::SSI1] = ssi;
//                    }
//                    else if(n == 2)
//                    {
//                        tvMap[TypeID::LLI2] = lli; tvMap[TypeID::SSI2] = ssi;
//                    }
//                    else if(n == 3)
//                    {
//                        tvMap[TypeID::LLI3] = lli; tvMap[TypeID::SSI3] = ssi;
//                    }
//                    else if(n == 5)
//                    {
//                        tvMap[TypeID::LLI5] = lli; tvMap[TypeID::SSI5] = ssi;
//                    }
//                    else if(n == 6)
//                    {
//                        tvMap[TypeID::LLI6] = lli; tvMap[TypeID::SSI6] = ssi;
//                    }
//                    else if(n == 7)
//                    {
//                        tvMap[TypeID::LLI7] = lli; tvMap[TypeID::SSI7] = ssi;
//                    }
//                    else if(n == 8)
//                    {
//                        tvMap[TypeID::LLI8] = lli; tvMap[TypeID::SSI8] = ssi;
//                    }
//                    else if(n == 9)
//                    {
//                        tvMap[TypeID::LLI9] = lli; tvMap[TypeID::SSI9] = ssi;
//                    }
//                }
//                else
//                {
//                    tvMap[type_2c] = data;
//                }
//            }
//
////           cout << endl;
//
//            theMap[sat] = tvMap;
//
//        }   // End loop over all the satellite
//
//        return theMap;
//    }

}  // End of namespace gpstk
