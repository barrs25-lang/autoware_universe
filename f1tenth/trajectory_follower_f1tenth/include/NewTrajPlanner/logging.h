#ifndef LOGGER_H
#define LOGGER_H

// Standard header files
#include <read_files.h>
#include <iomanip>
#include <chrono>
#include <time.h>

// This class is used to log data 
// for debugging or analysis purposes
class LOGGER {

	public:

		LOGGER();
		~LOGGER();

		// 1 pointer to a character implies you wish to write down a system message
		int writeToLog(const string& message);
		int writeToLog(const string& message, bool end_line);
		int openLogFile(const string& file_name);

	private:

		bool message_of_status = 0;
		ofstream message_of;

};

#endif