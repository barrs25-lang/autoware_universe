// Simulation environment for the DARPA Tactical Mapping Project
// Author: Julius Allen Marshall
// Date Created: December 6th, 2021
// Last Modified: December 6th, 2021
// Contact: mjulius@vt.edu

// File Decscription ####################################################################################
// This source file starts the simulated environment.
// End File Decscription ################################################################################

// Custom header files
#include "logging.h"

auto logger_time_start = std::chrono::high_resolution_clock::now();

// LOGGER constructor opens all log files. If a log fails to open, the program will reply to any call 
// to write to that log with an error.
LOGGER::LOGGER()
{


};

LOGGER::~LOGGER()
{

	if (message_of.is_open())
	{
		message_of.close();
	}
	
};

int LOGGER::openLogFile(const string& file_name)
{

  	message_of.open("../../Diagnostic_Logs/Dgnstc_" + file_name + "_Log.txt");

	if (!message_of) 
	{
		cout << "Unable to open file: ../../Diagnostic_Logs/Dgnstc_" + file_name + "_Log.txt" << endl;
		message_of_status = 0;
	}
	else
	{
		message_of_status = 1;
	}

}

// 1 pointer to a character and 1 double implies you wish to write down a system message with a timestamp
int LOGGER::writeToLog(const string& message, bool end_line)
{

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeElapsed = end_time - logger_time_start;

	int status = 0;

	if (message_of_status)
	{
		if (end_line)
		{
			message_of << message << std::endl;
		}
		else
		{
			message_of << message;
		}
		status = 1;
	}
	else
	{
		status = 0;
	}

	return status;

}

// 1 pointer to a character and 1 double implies you wish to write down a system message with a timestamp
int LOGGER::writeToLog(const string& message)
{

	auto end_time = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeElapsed = end_time - logger_time_start;

	int status = 0;

	if (message_of_status)
	{
		message_of << "TIME: " << setprecision(4) << timeElapsed.count();
		message_of << " - SYSTEM: " << message << std::endl;

		status = 1;
	}
	else
	{
		status = 0;
	}

	return status;

}