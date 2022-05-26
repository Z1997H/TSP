package cn.edu.jxufe.agl;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Random;
import java.util.Scanner;

public class TabuSearchTSP {
	
	public static final int MAX_GEN=1000;//���ĵ�������
	public static final int N=200;//ÿ����������ĸ���
	public static final int ll=20;//���ɳ���
	public static int cityNum=48;//�����������ֶ�����
	public static int jinji[][]=new int[ll][cityNum];//���ɱ�
	public static int[]Ghh=new int[cityNum];//��ʼ·������
	public static int bestGh[]=new int[cityNum];//��õ�·������
	public static int[]LocalGh=new int[cityNum];//��ǰ���·������
	public static int[]tempGh=new int[cityNum];//�����ʱ����
	public static int bestT;//��ѵĵ�������
	public static int bestEvaluation;
	public static int LocalEvaluation;
	public static int tempEvaluation;
	public static int point[][]=new int[cityNum][2];//ÿ�����е�����
	public static int dist[][]=new int[cityNum][cityNum];//�������
	public static int t;//��ǰ����
	public static Random random;
	
	//��ȡ���ݲ���ʼ��
	public static void read_data(String filepath) throws FileNotFoundException {
		String line=null;
		String substr[]=null;
		Scanner cin=new Scanner(new BufferedReader(new FileReader(filepath)));
		for (int i = 0; i < cityNum; i++) {
			line=cin.nextLine();
			line.trim();
			substr=line.split(" ");
			 point[i][0]=Integer.parseInt(substr[1]);//x����
			 point[i][1]=Integer.parseInt(substr[2]);//y����
		}
		cin.close();
		//����������ע������ļ��㷽ʽ�����õ���αŷʽ����
		for (int i = 0; i < cityNum; i++) {
			dist[i][i]=0;//�Խ���Ԫ��Ϊ0
			for (int j = i+1; j < cityNum; j++) {
				double rij=Math.sqrt((Math.pow(point[i][0]-point[j][0], 2)+
						             Math.pow(point[i][1]-point[j][1], 2))/10.0);
				//rij��������ȡ��
				int tij=(int) Math.round(rij);
				if(tij<rij) {
					dist[i][j]=tij+1;
					dist[j][i]=dist[i][j];
				}else {
					dist[i][j]=tij;
					dist[j][i]=dist[i][j];
				}
			}
		}
		dist[cityNum-1][cityNum-1]=0;
		
		t=0;
		bestT=0;
		bestEvaluation=Integer.MAX_VALUE;
		LocalEvaluation=Integer.MAX_VALUE;
		tempEvaluation=Integer.MAX_VALUE;
		random=new Random(System.currentTimeMillis());
	}
	
	
	//��ʼ������Ghh
	public static void initGroup() {
		Ghh[0]=random.nextInt(65535)%cityNum;
		int i,j;
		//ʹ�ò�����ÿ�����򶼲�һ��
		for (i = 1; i < cityNum;) {
			Ghh[i]=random.nextInt(65535)%cityNum;
			for (j = 0; j < i; j++) {
				if(Ghh[i]==Ghh[j]) {
					break;
				}
			}
			if(j==i) {
				i++;
			}
		}
	}
	
	//���Ʋ�������Gha���Ƶ�Ghb
	public static void copyGh(int[]Gha,int[]Ghb) {
		for (int i = 0; i < cityNum; i++) {
			Ghb[i]=Gha[i];
		}
	}
	
	//���ۺ���
	public static int evaluate(int[] chr) {
		int len=0;
		for (int i = 1; i < cityNum; i++) {
			len+=dist[chr[i-1]][chr[i]];
		}
		len+=dist[chr[cityNum-1]][chr[0]];
		return len;
	}
	
	//���򽻻����õ��ھ�
	public static void Linju(int[]Gh,int[]tempGh) {
		int i,temp;
		int rand1,rand2;
		//��Gh���Ƶ�tempGh
		for (i = 0; i < cityNum; i++) {
			tempGh[i]=Gh[i];
		}
		rand1=random.nextInt(65535)%cityNum;
		rand2=random.nextInt(65535)%cityNum;
		while(rand1==rand2) {
			rand2=random.nextInt(65535)%cityNum;
		}
		//����
		temp=tempGh[rand1];
		tempGh[rand1]=tempGh[rand2];
		tempGh[rand2]=temp;
	}
	
	//�жϱ����Ƿ��ڽ��ɱ���
	public static int panduan(int[]tempGh) {
		int i,j;
		int flag=0;
		for (i = 0; i <ll; i++) {
			flag=0;
			for ( j = 0; j < cityNum; j++) {
				if(tempGh[j]!=jinji[i][j]) {
					flag=1;//����ͬ
					break;
				}
			}
			if(flag==0) {	//��ͬ�����ش�����ͬ
			break;
			}
		}
		if(i==ll) {
			return 0;//������
		}else {
			return 1;//����
		}
	}
	
	//�������������
	public static void jiejinji(int[]tempGh) {
		int i,j,k;
		//ɾ�����ɱ��һ�����룬���������ǰ�ƶ�
		for (i = 0; i < ll-1; i++) {
			for (j = 0; j < cityNum; j++) {
				jinji[i][j]=jinji[i+1][j];
			}
		}
		//�µı��������ɱ�
		for (k = 0; k < cityNum; k++) {
			jinji[ll-1][k]=tempGh[k];
		}
	}
	
	public static void solve() {
		int nn;
		//��ʼ������Ghh
		initGroup();
		copyGh(Ghh,bestGh);//���Ƶ�ǰ����Ghh����ñ���bestGh
		bestEvaluation=evaluate(Ghh);
		
		while(t<MAX_GEN) {
			nn=0;
			LocalEvaluation=Integer.MAX_VALUE;
			while(nn<N) {
				Linju(Ghh,tempGh);//�õ���ǰ����Ghh���ھӱ���tempGh
				if(panduan(tempGh)==0) {//�ж��Ƿ��ڽ��ɱ���
					//����
					tempEvaluation=evaluate(tempGh);
					if(tempEvaluation<LocalEvaluation) {
						copyGh(tempGh,LocalGh);
						LocalEvaluation=tempEvaluation;
					}
					nn++;
				}
			}
			if(LocalEvaluation<bestEvaluation) {
				bestT=t;
				copyGh(LocalGh,bestGh);
				bestEvaluation=LocalEvaluation;
			}
			copyGh(LocalGh,Ghh);
			
			//����ɱ�LocalGh������ɱ�
			jiejinji(LocalGh);
			t++;
		}
		System.out.println("��ѵ�������:"+bestT);
		System.out.println("��ѳ���Ϊ:"+bestEvaluation);
		System.out.println("���·��Ϊ:");
		for (int i = 0; i < cityNum; i++) {
			System.out.print(bestGh[i]+"-->");
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException {
		System.out.println("Start....");
		long start = System.currentTimeMillis();
		TabuSearchTSP.read_data("./src/cn/edu/jxufe/agl/data.txt");
		TabuSearchTSP.solve();
		long end = System.currentTimeMillis();
		System.out.println("\n��ʱ��" + (end-start) + "ms");
	}
}