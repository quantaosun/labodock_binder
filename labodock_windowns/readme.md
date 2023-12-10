## Both Windows PowerShell and Docker-Desktop can do the docking in our case

### **But if you feel not comfortable with terminals, here is a quick tutorial for Docker-Desktop**

### Step 1: Search and pull the docking image from the Docker hub

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/45a8c3c7-ebd2-40aa-b206-c6afca5d69cc)

### Step 2: Click the play button of the pulled docking image

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/6f1718ff-6561-412b-a62e-081fa5805471)

### Step 3: Do not click **Run** directly, but open the dropdown menu

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/89f48550-e862-4865-85e2-8848feb540da)

### Step 4: Fill the port as ```8888```, and fill the **container path** as ```/home/jovyan/work```, these two should remain unchanged. Copy the folder path where you want to store the docking result and paste it inside **volume path**

This step is critical, Double-check that everything is good, and pay attention to the two paths, The volume path is a Windows path of your local folder, and the container path is a Linux path that comes with a different style of format. /home/jovyan/work, the name jovyan is just a jupyter notebook community member name, just keep it like that, since no matter who starts this docker container, the path of this container is always named like this.

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/a80073f9-51eb-4346-af8d-0817c27a7d77)

### Step 5: Click Run and then click the 8888:8888 link on the top of the interface and your browser should open now

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/3163ed09-31ab-4da2-966c-f6367d6b1d41)

### Step 6: The browser opens a jupyter notebook and asks for a token pin

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/29492fd0-856c-4d7d-ba57-7e52c18df597)

### Step 7: Go back to the Docker 8888:8888 interface, and find the token string beneath, copy exactly how long it is

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/cf0289af-2b34-4a66-aae4-689ee2647373)

### Step 8: Go to the browser/Jupyter interface and paste this token string and click login

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/b3d83de2-8c42-482b-81a1-dfc3cac86201)

### Step 9: Double-click ```wedock_win.ipynb``` and once it opens, you are ready to do docking on Windows locally now

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/8eb29b2f-3dd6-4ad6-9fdc-7f25137699b2)

### Step 10: Start docking, after having modified the three variables required, you click **Run** > **Run all**

(The top right kernel should always be selected as **notebook** all the time, normally it will do it by default)

![image](https://github.com/quantaosun/labodock_binder/assets/75652473/a297a000-9fb2-4d18-bd0d-249bdbfd7743)








