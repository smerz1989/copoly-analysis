import os
import paramiko as pk
import re
import sys
from pssh.clients import ParallelSSHClient as psshclient
from pssh.clients import SSHClient as psshclient
from gevent import joinall

class ServerConnection(object):
    def __init__(self):
        self.ssh_client = self.connect_to_server()

    def connect_to_server(self):
        #ssh_client = pk.SSHClient()
        #ssh_client.set_missing_host_key_policy(pk.AutoAddPolicy())
        try:
            username = os.environ['COPOLYUSERNAME']
            password = os.environ['COPOLYPASSWORD']
            hostname = os.environ['COPOLYHOST']
            self.host = hostname
        except KeyError:
            print('COPOLYHOST, COPOLYUSERNAME, or COPOLYPASSWORD environment variable is not set.  Set these environment variables to the correct values for target server to fix issues.')
        #ssh_client.connect(hostname=hostname,username=username,password=password)
        ssh_client = psshclient(hostname,user=username,password=password)
        return ssh_client

    def exec_command(self,command):
        output = self.ssh_client.run_command(command)
        print(output)
        return((output[4],output[2],output[3]))

    def check_if_file_exists(self,remote_path):
        stdin, stdout, stderr = self.exec_command('ls '+str(remote_path))
        #stdin, stdout, stderr = self.ssh_client.exec_command('ls '+str(remote_path))
        print(list(stdout))
        directory_check = '\n'.join(list(stderr))
        if len(re.findall(r'No\ such\ file',directory_check))>0:
            return False
        else:
            return True
            
    def mkdir(self,path):
        #Check if directory already exists
        if self.check_if_file_exists(path): 
            print("Directory already exists")
            raise OSError
        else:
            self.exec_command('mkdir -p '+str(path))
            #self.ssh_client.exec_command('mkdir -p '+str(path))

    def lsdir(self,path):
        stdin,stdout,stderr = self.exec_command('ls '+path)  
        #ftp_client = self.ssh_client.open_sftp()
        return(list(stdout))

    def send_file(self,local_path,remote_path):
        #ftp_client = self.ssh_client.open_sftp()
        print("Transferring file: {} to server".format(local_path))
        self.ssh_client.scp_send(local_path,remote_path)
        #greenlets = self.ssh_client.copy_file(local_path,remote_path)
        #joinall(greenlets,raise_error=True) 
        #ftp_client.put(local_path,remote_path,callback=self.print_progress)
        #ftp_client.close()

    def get_file(self,remote_path,local_path):
        #ftp_client = self.ssh_client.open_sftp()
        print("Downloading file: {} from server".format(local_path))
        self.ssh_client.scp_rcv(remote_path,local_path)
        #greenlets = self.ssh_client.copy_remote_file(remote_path,local_path)
        #joinall(greenlets,raise_error=True)
        #ftp_client.get(remote_path,local_path,callback=self.print_progress)
        #sys.stdout.write('\n')
        #sys.stdout.flush()
        #ftp_client.close()

    def print_progress(self,transferred, toBeTransferred):
        progress = transferred/int(toBeTransferred)*100
        #print("{}% transferred of {} bytes to be transferred".format(progress,toBeTransferred))
        sys.stdout.write('\rTransfer Progress: %.2f %%' % (progress))
        sys.stdout.flush()
