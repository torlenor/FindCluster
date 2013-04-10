#ifndef THREEDCLUSTERS_KEYS_HPP
#define THREEDCLUSTERS_KEYS_HPP

void setHightlightCluster(int cnt){
	red=sred;
	green=sgreen;
	blue=sblue;
	for(int i1=0;i1<leng1;i1++)
	for(int i2=0;i2<leng2;i2++)
	for(int i3=0;i3<leng3;i3++){
		int is = i1 + i2*leng1 + i3*leng1*leng2;
		if(lpoints.at(i1).at(i2).at(i3)==cnt){
			red[is]=1; green[is]=1; blue[is]=1;
		}
	}
}

void setOnlyCluster(int cnt){
	for(int i1=0;i1<leng1;i1++)
	for(int i2=0;i2<leng2;i2++)
	for(int i3=0;i3<leng3;i3++){
		int is = i1 + i2*leng1 + i3*leng1*leng2;
		if(lpoints.at(i1).at(i2).at(i3)!=cnt){
		red[is]=0; green[is]=0; blue[is]=0;
		}else{
			red[is]=1; green[is]=1; blue[is]=1;
		}
	}
	
}

void setSectorColors(bool set){
	if(set==true){
		for(int i1=0;i1<leng1;i1++)
		for(int i2=0;i2<leng2;i2++)
		for(int i3=0;i3<leng3;i3++){
			int is = i1 + i2*leng1 + i3*leng1*leng2;
			if(isinsector[is] == -1){
				red[is]=1; green[is]=0; blue[is]=0;
			}else if(isinsector[is] == 0){
				red[is]=0; green[is]=1; blue[is]=0;
			}else if(isinsector[is] == 1){
				red[is]=0; green[is]=0; blue[is]=1;
			}
		}
	}else{
		red=sred;
		green=sgreen;
		blue=sblue;
	}
}

void setOnlySector(int sec){
	for(int i1=0;i1<leng1;i1++)
	for(int i2=0;i2<leng2;i2++)
	for(int i3=0;i3<leng3;i3++){
		int is = i1 + i2*leng1 + i3*leng1*leng2;
		if(isinsector[is]!=sec){
			red[is]=0; green[is]=0; blue[is]=0;
		}else{
			if(sec==-1){
				red[is]=1; green[is]=0; blue[is]=0;
			}else if(sec==0){
				red[is]=0; green[is]=1; blue[is]=0;
			}else if(sec==1){
				red[is]=0; green[is]=0; blue[is]=1;
			}
		}
	}
}

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
			setHightlightCluster(cnt);
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			onecluster=0;
			cnt++;
			if(cnt > nclusters-1)
				cnt = 0;
			setHightlightCluster(cnt);
			cout << "Cluster " << cnt << " selected!" << endl;
		}
	}else if(key == 's' || key == 'S' ){
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			onecluster=1;
			cnt--;
			if(cnt < 0)
				cnt = nclusters-1;
			setOnlyCluster(cnt);
			cout << "Cluster " << cnt << " selected!" << endl;
		}else{
			onecluster=1;
			cnt++;
			if(cnt > nclusters-1)
				cnt = 0;
			setOnlyCluster(cnt);
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
		setOnlySector(-1);
	}else if(key == '2'){
		// Select the sector 0
		setOnlySector(0);
	}else if(key == '3'){
		// Select the sector 1
		setOnlySector(1);
	}else if(key == 'z' || key == 'Z'){
		// Sector colors
		int mod = glutGetModifiers();
		if (mod == GLUT_ACTIVE_SHIFT){
			setSectorColors(false);
		}else{
			setSectorColors(true);
		}
	}
}

#endif //THREEDCLUSTERS_KEYS_HPP
