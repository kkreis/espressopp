#ifndef PyLogger_HPP
#define PyLogger_HPP

/** \file PyLogger.hpp    Python implementation for logging.

<b>Responsible:</b>
<a href="mailto:brandes@scai.fraunhofer.de">Thomas Brandes</a>

*/

#include <esutil/Logger.hpp>
#include "python.hpp"

namespace log4espp {

  /**************************************************************************
  *                                                                         *
  *                                                                         *
  *                                                                         *
  **************************************************************************/

  /** PyLogger is an implementaton of the abstract class Logger that uses 
      corresponding Python Logger objects for the output of the logging 
      statements.

      One major efficiency aspect is that an object of this class will not
      ask for the logging level of its Python equivalent; therefore it has 
      to be ensured that changing the logging level of a Python logger
      will notity an object of this C++ class.

  */

  class PyLogger : public Logger {

  private:

   boost::python::object pyLogger;   //!< Pointer to the Python instance of this logger.

   /** This routine updates the logging level by the given level of the Python logger

       \param pyLevel is the Python coding of the level

       The Python coding of the level will be translated to Logger::Level.
   */

   void setPythonLevel(int pyLevel);

   /** Commonly used routine for logging output that can be used for all levels

       \param level is the string representation of the level
       \param loc is the file location of the logging statement
       \param msg is the message output of the logging statement

       The logging output will be printed via the Python Logger object. 
   */

   void log(int level, Location& loc, const std::string& msg);

  public:

   typedef boost::python::object object;

   /** Initialization routine to be called with initialization */

   static void initLogging();

   static void registerPython();

   /** Constructor for a Logger that will use a corresponding Python instance
       for the configuration and for the output of the logging messages  
   */

   PyLogger(std::string, class Logger* parent);

   ~PyLogger() {}

   /** This routine sets the Python instance of the logger and/or 
       updates the logging level of the C++ class.

       \param pyLogger is the Python Logger object.

       This routine must also be called to update the logging level
       of this object.
   */

   void setPythonLogger(boost::python::object pyLogger);

   /** This routine sets the Python loggers for this logger and the
       descendants.

       \param ParentName is the name of the parent of this logger.

   */

   void setPythonLoggers(std::string& parentName);

   /** Implementation of Logger::trace */

   virtual void trace(Location loc, const std::string& msg);

   /** Implementation of Logger::debug */

   virtual void debug(Location loc, const std::string& msg);

   /** Implementation of Logger::info */

   virtual void info (Location loc, const std::string& msg);

   /** Implementation of Logger::warn */

   virtual void warn (Location loc, const std::string& msg);

   /** Implementation of Logger::error */

   virtual void error(Location loc, const std::string& msg);

   /** Implementation of Logger::fatal */

   virtual void fatal(Location loc, const std::string& msg);

  };

}

#endif
