# Rockfall-Simulator

**Basic Information**

The directory provides the raw code (.m and .fig files) which need to be compiled in MatLab in order to become a useful program.

Included in the directory is the actual .exe file (Falling.exe) which when executed opens the program GUI. 
Note that to open the program you must have MATLAB Compiler Runtime (MCR) installed running version 8.0 (R2012b)



**Program Description**

This is a very basic rockfall simulator done in MatLab as a school project. 

The simulator assumes rocks to be point objects and launches them down a pre-defined slope. What has been coded is the projectile algorithm, impact algorithm, and rolling algorithm.

Slope parameters and other input properties are entered into a very basic GUI for the code to execute. 

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/7b3bcf9e-596f-4b0f-9eb9-1c838f8362c9)

Prior knowledge of rockfall programs and their input parameters is required to oeprate the program. The program follows a basic decision flowchart as shown below.

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/107fb5bd-2a20-4c03-8436-8e1a84174260)




**Program Limitations**

Use of this program is entirely at your own risk and there is no warranty whatsoever provided.

There is a limitation in the code of the program such that rocks that start rolling cannot enter back into the projectile and impact algorithm. The user must manually gather the velocity and other data just at the point where the rocks would launch off the cliff (shown below) and re-start the code. 

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/87f1cf30-6715-437a-b0b8-c44134b04434)

Once the vector components of velocity are collected as shown below

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/119a9502-f74e-4714-bd5f-fb41c3023e81)
![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/96f5e0bb-fad7-44aa-87d1-d766b14088de)

The program is re-executed and proper results are obtained.

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/d43adbe6-5306-4130-8fd0-788f1c534e52)

**Verification**

Some basic verifications were done against RocFall (Rocscience program). Rocks of 1000kg mass & 0.446m radius were considered. Note that at the time of writing in 2014, RocFall point object algorithm considered rocks to "slide" downhill instead of "roll". Therefore kinetic friction and not rolling friction is used. Kinetic friction is derived from a friction angle formula as shown below & conversion was done.

Kinetic Friction = tan (friction angle)

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/2c62875c-ce5c-4bcd-a8c0-8a82fdd6887d)



**Verification 1**

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/0628f4a2-4432-4087-ae8a-a87d57cb0d7a)



**Verification 2**

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/dd520907-c89f-4b4e-98f1-f7ac7c6b3d4d)



**Verification 3**

![image](https://github.com/akatragjini/Rockfall-Simulator/assets/51655773/42e407d2-8750-4347-998b-957954832bdc)


