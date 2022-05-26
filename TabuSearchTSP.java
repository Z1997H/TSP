package cn.edu.jxufe.agl;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.Random;
import java.util.Scanner;

public class TabuSearchTSP {
	
	public static final int MAX_GEN=1000;//最大的迭代次数
	public static final int N=200;//每次搜索领域的个数
	public static final int ll=20;//禁忌长度
	public static int cityNum=48;//城市数量，手动设置
	public static int jinji[][]=new int[ll][cityNum];//禁忌表
	public static int[]Ghh=new int[cityNum];//初始路径编码
	public static int bestGh[]=new int[cityNum];//最好的路径编码
	public static int[]LocalGh=new int[cityNum];//当前最好路径编码
	public static int[]tempGh=new int[cityNum];//存放临时编码
	public static int bestT;//最佳的迭代次数
	public static int bestEvaluation;
	public static int LocalEvaluation;
	public static int tempEvaluation;
	public static int point[][]=new int[cityNum][2];//每个城市的坐标
	public static int dist[][]=new int[cityNum][cityNum];//距离矩阵
	public static int t;//当前迭代
	public static Random random;
	
	//读取数据并初始化
	public static void read_data(String filepath) throws FileNotFoundException {
		String line=null;
		String substr[]=null;
		Scanner cin=new Scanner(new BufferedReader(new FileReader(filepath)));
		for (int i = 0; i < cityNum; i++) {
			line=cin.nextLine();
			line.trim();
			substr=line.split(" ");
			 point[i][0]=Integer.parseInt(substr[1]);//x坐标
			 point[i][1]=Integer.parseInt(substr[2]);//y坐标
		}
		cin.close();
		//计算距离矩阵，注意这里的计算方式，才用的是伪欧式距离
		for (int i = 0; i < cityNum; i++) {
			dist[i][i]=0;//对角线元素为0
			for (int j = i+1; j < cityNum; j++) {
				double rij=Math.sqrt((Math.pow(point[i][0]-point[j][0], 2)+
						             Math.pow(point[i][1]-point[j][1], 2))/10.0);
				//rij四舍五入取整
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
	
	
	//初始化编码Ghh
	public static void initGroup() {
		Ghh[0]=random.nextInt(65535)%cityNum;
		int i,j;
		//使得产生的每个基因都不一样
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
	
	//复制操作，将Gha复制到Ghb
	public static void copyGh(int[]Gha,int[]Ghb) {
		for (int i = 0; i < cityNum; i++) {
			Ghb[i]=Gha[i];
		}
	}
	
	//评价函数
	public static int evaluate(int[] chr) {
		int len=0;
		for (int i = 1; i < cityNum; i++) {
			len+=dist[chr[i-1]][chr[i]];
		}
		len+=dist[chr[cityNum-1]][chr[0]];
		return len;
	}
	
	//领域交换，得到邻居
	public static void Linju(int[]Gh,int[]tempGh) {
		int i,temp;
		int rand1,rand2;
		//将Gh复制到tempGh
		for (i = 0; i < cityNum; i++) {
			tempGh[i]=Gh[i];
		}
		rand1=random.nextInt(65535)%cityNum;
		rand2=random.nextInt(65535)%cityNum;
		while(rand1==rand2) {
			rand2=random.nextInt(65535)%cityNum;
		}
		//交换
		temp=tempGh[rand1];
		tempGh[rand1]=tempGh[rand2];
		tempGh[rand2]=temp;
	}
	
	//判断编码是否在禁忌表中
	public static int panduan(int[]tempGh) {
		int i,j;
		int flag=0;
		for (i = 0; i <ll; i++) {
			flag=0;
			for ( j = 0; j < cityNum; j++) {
				if(tempGh[j]!=jinji[i][j]) {
					flag=1;//不相同
					break;
				}
			}
			if(flag==0) {	//相同，返回存在相同
			break;
			}
		}
		if(i==ll) {
			return 0;//不存在
		}else {
			return 1;//存在
		}
	}
	
	//解禁忌与加入禁忌
	public static void jiejinji(int[]tempGh) {
		int i,j,k;
		//删除禁忌表第一个编码，后面编码往前移动
		for (i = 0; i < ll-1; i++) {
			for (j = 0; j < cityNum; j++) {
				jinji[i][j]=jinji[i+1][j];
			}
		}
		//新的编码加入禁忌表
		for (k = 0; k < cityNum; k++) {
			jinji[ll-1][k]=tempGh[k];
		}
	}
	
	public static void solve() {
		int nn;
		//初始化编码Ghh
		initGroup();
		copyGh(Ghh,bestGh);//复制当前编码Ghh到最好编码bestGh
		bestEvaluation=evaluate(Ghh);
		
		while(t<MAX_GEN) {
			nn=0;
			LocalEvaluation=Integer.MAX_VALUE;
			while(nn<N) {
				Linju(Ghh,tempGh);//得到当前编码Ghh到邻居编码tempGh
				if(panduan(tempGh)==0) {//判断是否在禁忌表中
					//不在
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
			
			//解禁忌表，LocalGh加入禁忌表
			jiejinji(LocalGh);
			t++;
		}
		System.out.println("最佳迭代次数:"+bestT);
		System.out.println("最佳长度为:"+bestEvaluation);
		System.out.println("最佳路径为:");
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
		System.out.println("\n耗时：" + (end-start) + "ms");
	}
}