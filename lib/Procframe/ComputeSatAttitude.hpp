
#ifndef COMPUTE_SAT_ATTITUDE_HPP
#define COMPUTE_SAT_ATTITUDE_HPP

#include "ProcessingClass.hpp"
#include "ReferenceSystem.hpp"
#include "SolarSystem.hpp"
#include "AntexReader.hpp"


namespace gpstk
{

    struct YawData
    {
        YawData()
            : type("UNK"), beta(0.0),
              event(0), nominal(0.0), modeled(0.0)
        {};

        std::string type;
        double beta;
        int event;
        double nominal;
        double modeled;
    };

    typedef std::map<SatID,YawData> satYawDataMap;


    class ComputeSatAttitude
    {
    public:

        /// Default constructor.
        ComputeSatAttitude()
            : pRefSys(NULL), pSolSys(NULL), pAntexReader(NULL)
        {};


        /** Set ReferenceSystem object to be used.
        *
        * @param ref  ReferenceSystem object.
        */
        virtual ComputeSatAttitude& setReferenceSystem(ReferenceSystem& ref)
        { pRefSys = &ref; return (*this); };


        /** Set SolarSystem object to be used.
        *
        * @param sol  SolarSystem object.
        */
        virtual ComputeSatAttitude& setSolarSystem(SolarSystem& sol)
        { pSolSys = &sol; return (*this); };


        /** Set AntexReader object to be used.
        *
        * @param antex  AntexReader object.
        */
        virtual ComputeSatAttitude& setAntexReader(AntexReader& antex)
        { pAntexReader = &antex; return (*this); };


        satYawDataMap getAttitude( const CommonTime& epoch,
                                   const satVectorMap& orbits )
            throw(ProcessingException);


        /// Return a string identifying this object.
        virtual std::string getClassName(void) const;


        /// Destructor
        virtual ~ComputeSatAttitude() {};

    private:

        bool noonManeuver(const CommonTime& epoch,
                          const double beta,
                          const double betaLimit,
                          const double mu,
                          const double muRate,
                          const double maxYawRate,
                          const double yawBias,
                          const Vector<double>& unit_rsat,
                          const Vector<double>& unit_rsun,
                          double& modeledYaw);

        bool nightManeuver(const CommonTime& epoch,
                           const double beta,
                           const double betaLimit,
                           const double mu,
                           const double muRate,
                           const double maxYawRate,
                           const double yawBias,
                           const Vector<double>& unit_rsat,
                           const Vector<double>& unit_rsun,
                           double& modeledYaw);

        bool shadowManeuver(const CommonTime& epoch,
                            const double beta,
                            const double betaLimit,
                            const double mu,
                            const double muRate,
                            const double maxYawRate,
                            const double yawBias,
                            const Vector<double>& unit_rsat,
                            const Vector<double>& unit_rsun,
                            double& modeledYaw);

    private:

        /// Pointer to ReferenceSystem
        ReferenceSystem* pRefSys;

        /// Pointer to SolarSystem
        SolarSystem* pSolSys;

        /// Pointer to AntexReader
        AntexReader* pAntexReader;

        double nominalYaw;
        double modeledYaw;

        satYawDataMap satYawData;

    }; // End of class 'ComputeSatAttitude'

    //@}

    }  // End of namespace gpstk

#endif   // COMPUTE_SAT_ATTITUDE_HPP
