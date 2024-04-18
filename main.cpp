/*********************************************
Plane Frame Analysis - Main program
*********************************************/
#include "frame.h"
#include "clockEXH.h"

int main (int argc, char *argv[])
{
    // CArrayBase class is used to track memory allocation, deallocation
    // and route the error messages from the CVector and CMatrix classes
    CArrayBase AB;
    {
        CFrame TheFrame; // the one and only frame!
        try
        {
	        // show program banner
            TheFrame.Banner (std::cout);

	        // Prepare for I/O
	        TheFrame.PrepareIO (argc, argv);
	
            // start the timer --------------------------------------------------------
            CClock Timer;
            std::string strDateTime;
            Timer.GetDateTime (strDateTime);
            std::cout << "\nStarting out at : " << strDateTime << "\n";

            // read the data and analyze
            TheFrame.Analyze ();

            // end the timer --------------------------------------------------------
            // get the current date and time
            Timer.GetDateTime (strDateTime);
            std::cout << "\n      Ending at : " << strDateTime << "\n";
            // compute the elapsed time -----------------------------------------------
            std::cout << "Elapsed wall clock time: " << Timer.DiffTime ()
                      << " seconds\n";
        }
        // --------------------
        // trap all errors here
        // --------------------
        // errors trapped by this program
        catch (CLocalErrorHandler::ERRORCODE &err)
        {
            TheFrame.DisplayErrorMessage (err);
        }

        // errors trapped by the library functions
        catch (CGlobalErrorHandler::ERRORCODE &err)
        {
            CGlobalErrorHandler::ErrorHandler (err);
        }

        // errors trapped in the CVector/CMatrix classes
        catch (CArrayBase::ERRORVM &err)
        {
            CGlobalErrorHandler::ErrorHandler (err);
        }

        // errors trapped by C++ 
        catch (std::exception &err)
        {
            CGlobalErrorHandler::ErrorHandler (err);
        }

        // forgotten to handle a trapped error?
        catch (...)
        {
            std::cout << "Sorry, could not catch the error whatever it is.\n";
        }

        // Close input and output files
        TheFrame.TerminateProgram ();
    }

    AB.ShowStatistics (std::cout);

	return 0;
}