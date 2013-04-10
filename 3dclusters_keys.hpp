#ifndef THREEDCLUSTERS_KEYS_HPP
#define THREEDCLUSTERS_KEYS_HPP

void processNormalKeys(unsigned char key, int x, int y){
	if(key == 27){
		exit(0);
	}else if(key == '+'){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_ALT){
			 alpha += 0.1;
			 if(alpha > 1.0)
			 	alpha=1.0;
		}
		else{
			sphereradius += 0.1;
			pointsize += 1.0;
			if(sphereradius > 1.0)
				sphereradius = 1.0;
		}
		cout << "Radius = " << pointsize << " Alpha = "  << alpha << endl;
	}else if(key == '-'){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_ALT){
			 alpha -= 0.1;
			 if(alpha < 0.0)
			 	alpha=0;
		}
		else{
			sphereradius -= 0.1;
			pointsize -= 1.0;
			if(pointsize < 0.0)
				pointsize = 1.0;
			if(sphereradius < 0.0)
				sphereradius = 0.0;
		}
		cout << "Radius = " << pointsize << " Alpha = "  << alpha << endl;
	}else if(key == 'c' || key == 'C' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			onecluster=0;
			cnt--;
			if(cnt < 0)
				cnt = nclusters-1;
			red=sred;
			green=sgreen;
			blue=sblue;
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			onecluster=0;
			cnt++;
			if(cnt > nclusters-1)
				cnt = 0;
			red=sred;
			green=sgreen;
			blue=sblue;
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}
	}else if(key == 's' || key == 'S' ){

		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			onecluster=1;
			cnt--;
			if(cnt < 0)
				cnt = nclusters-1;
			for(int i=0;i<nclusters;i++){
				red[i]=0;
				green[i]=0;
				blue[i]=0;
			}
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			onecluster=1;
			cnt++;
			if(cnt > nclusters-1)
				cnt = 0;
			for(int i=0;i<nclusters;i++){
				red[i]=0;
				green[i]=0;
				blue[i]=0;
			}
			red[cnt]=1; green[cnt]=1; blue[cnt]=1;
			cout << "Cluster " << cnt << " selected!" << endl;
		}	
	}else if(key == 'r'){
		cnt=-1;
		red=sred;
		green=sgreen;
		blue=sblue;
		cout << "No Cluster selected!" << endl;
	}else if(key == 'n' || key == 'N' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			selconfig--;
			if(selconfig==-1)
				selconfig=nconfig-1;
			cout << "Loading configuration " << selconfig << " ..." << endl;
			cluster3input(selconfig);
		}else{
			selconfig++;
			if(selconfig==nconfig)
				selconfig=0;
			cout << "Loading configuration " << selconfig << " ..." << endl;
			cluster3input(selconfig);
		}
	}else if(key == '1'){
		// Select the sector -1

	}else if(key == '2'){
		// Select the sector 0
	
	}else if(key == '3'){
		// Select the sector 1
	
	}else if(key == 'c'){
		// Sector colors
		
	}
}

#endif //THREEDCLUSTERS_KEYS_HPP
