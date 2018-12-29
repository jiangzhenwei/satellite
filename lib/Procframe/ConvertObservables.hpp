#pragma ident "$ID$"

#ifndef CONVERT_OBSERVABLES_HPP
#define CONVERT_OBSERVABLES_HPP

#include "ProcessingClass.hpp"

namespace gpstk
{

    class ConvertObservables : public ProcessingClass
    {
    public :

        // Default constructr
        ConvertObservables() {};


        /** Return a satTypeValueMap object, adding the new data generated
         *  when calling this object.
         *
         * @param  time     Epoch corresponding to the data.
         * @param  gData    Data object holding the data.
         */
        virtual satTypeValueMap& Process( const CommonTime& time,
                                          satTypeValueMap& gData )
            throw(ProcessingException);


        /** Return a gnssSatTypeValue object, adding the new data
         *  generated when calling this object.
         *
         * @param gData     Data object holding the data.
         */
        virtual gnssSatTypeValue& Process(gnssSatTypeValue& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


        /** Return a gnssRinex object, adding the new data generated
         *  when calling this object.
         *
         * @param gData     Data object holding the data.
         */
        virtual gnssRinex& Process(gnssRinex& gData)
            throw(ProcessingException)
        { Process(gData.header.epoch, gData.body); return gData; };


        /** Return a gnssDataMap object, adding the new data generated
         *  when calling this object.
         *
         * @param gData     Data object holding the data.
         */
        virtual gnssDataMap& Process(gnssDataMap& gData)
            throw(ProcessingException);


        /** Method to set the satSys.
         *
         * @param sys      SatSys
         */
        virtual ConvertObservables& setSatSystem(const SatID::SatelliteSystem& sys)
        { satSys = sys; return (*this); };


        /** Method to get the satSys.
         */
        virtual SatID::SatelliteSystem getSatSystem() const
        { return satSys; };


        /** Method to set the oldType.
         *
         * @param type      Type
         */
        virtual ConvertObservables& setOldType(const TypeID& type)
        { oldType = type; return (*this); };


        /** Method to get the oldType.
         */
        virtual TypeID getOldType() const
        { return oldType; };


        /** Method to set the newType.
         *
         * @param type      Type
         */
        virtual ConvertObservables& setNewType(const TypeID& type)
        { newType = type; return (*this); };


        /** Method to get the newType.
         */
        virtual TypeID getNewType() const
        { return newType; };


        // Return a string identifying this object.
        virtual std::string getClassName () const;


        // Default deconstructr
        ~ConvertObservables() {};


    private:

        SatID::SatelliteSystem satSys;

        TypeID oldType;

        TypeID newType;

    }; // End of class 'ConvertObservables'

    // @}

}; // end of namespace gpstk

#endif   // CONVERT_OBSERVABLES_HPP
