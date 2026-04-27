#ifndef DATASEND_H_
#define DATASEND_H_

#include <fcntl.h>   
#include <unistd.h>
#include <string.h>
#include <sys/socket.h>	//unix socket
#include <sys/types.h>
#include <stdlib.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <sched.h>
#include <termios.h> 

class DataSend {
public:

	int sockfd = 0, Clientsockfd = 0, err;
	struct sockaddr_in serverInfo,clientInfo;
	unsigned int addrlen = sizeof(clientInfo);
	int ret = -1;
	char message[30] = "Hi,this is server.\n";
	// Socket declared
	DataSend() 
	{
		printf("Socket called!\n");
	}
	//check if there some error
	void init(uint16_t port_num) 
	{

		sockfd = socket(AF_INET, SOCK_STREAM, 0);
		if(sockfd == -1){
			printf("Fail to create a socket!\n");
		}
		//************************
		bzero(&serverInfo,sizeof(serverInfo));
		serverInfo.sin_family = AF_INET;
		// local host test
		serverInfo.sin_addr.s_addr = inet_addr("127.0.0.1");	//127.0.0.1		104.39.160.46	104.39.89.30
		serverInfo.sin_port = htons(port_num);

		std::cout << "Waiting to connect ..." << std::endl;
		while(1)
		{
			err = connect(sockfd,(struct sockaddr *)&serverInfo,sizeof(serverInfo));
	    	if(err==-1)
	    	{

			}
			else
			{
				printf("Connected!\n");
				std::cout << ">> Port:" << port_num << std::endl;
				break;
			}
		}
	}

	int set_accept()
	{
		printf("wait for connection....\n");
		Clientsockfd = accept(sockfd,(struct sockaddr*)&serverInfo, &addrlen);
		if (Clientsockfd == -1)
		{
			perror("Count not accept gnc connection");
		}
		else
		{
			printf("Connection Accepted!\n");
		}
		return Clientsockfd;
	}
	void set_nonblocking()
	{
		//fcntl(Clientsockfd , F_SETFL, O_NONBLOCK);
		int flag = fcntl(sockfd, F_GETFL, 0);
		if (flag < 0) {
			perror("fcntl1 F_GETFL fail");
			//return 0;
		}
		if (fcntl(sockfd, F_SETFL, flag | O_NONBLOCK) < 0) {
			perror("fcntl1 F_SETFL fail");
			//return 0;
		}
			printf("Non-blocking Activated!\n");
	}
	void process_sending(void * message, long unsigned int SIZE) 
	{

		//std::cout << "char size:" << sizeof(data) << std::endl;
		// if (flag == 0)
		// {
			// send(sockfd, message, SIZE, MSG_WAITALL);
			
		// }
		// else
		// {
			send(sockfd, message, SIZE, MSG_DONTWAIT);

		// }
		//send(sockfd, message, SIZE, MSG_DONTWAIT);
		// force cast
		//float(&pArray)[n_update][4] = *reinterpret_cast<float(*)[n_update][4]>(data);

	}
	// data receive from server
	int process_receiving(void * recvbuf, int SIZE, int flag) 
	{

		if (flag == 0)
		{
			ret = recv(sockfd, recvbuf, SIZE, MSG_WAITALL);
		}
		else if (flag == 1)
		{
			ret = recv(sockfd, recvbuf, SIZE, MSG_DONTWAIT);
		}
		return ret;

	}
	void End() 
	{

		if (err != -1)
		{
			shutdown(sockfd, SHUT_RDWR);
			close(sockfd);
		}
		printf("Socket Closed!\n");
	}
};

class DataSendWrapper
{

public: 
	DataSend* _in_pose_sim_socket;
	DataSend* _in_map_sim_socket;
	DataSend* _in_goal_generation_socket;
	DataSend* _out_path_sim_socket;
	DataSend* _out_path_mpc_socket;

	void socketWrapper(DataSend* _in_pose_sim_socket, DataSend* _in_map_sim_socket, DataSend* _in_goal_generation_socket, DataSend* _out_path_sim_socket, DataSend* _out_path_mpc_socket)
	{

		this->_in_pose_sim_socket = _in_pose_sim_socket;
		this->_in_map_sim_socket = _in_map_sim_socket;
		this->_in_goal_generation_socket = _in_goal_generation_socket;
		this->_out_path_sim_socket = _out_path_sim_socket;
		this->_out_path_mpc_socket = _out_path_mpc_socket;

	}

	DataSend* operator[](int i)
	{

		switch(i)
		{

			case 1: return _in_pose_sim_socket; break;
			case 2: return _in_map_sim_socket; break;
			case 3: return _in_goal_generation_socket; break;
			case 4: return _out_path_sim_socket; break;
			case 5: return _out_path_mpc_socket; break;

		}

	}

private:

};
	

#endif