void updateWeights(pcl::PointCloud<PointT>::Ptr &map, Database query, double kde_singma[], int t) {
	vector<double> dis_vec;
	double std_x = kde_singma[0];
	double std_y = kde_singma[1];
	pcl::search::KdTree<PointT>::Ptr kdtree(new pcl::search::KdTree<PointT>);
	// 利用KDTree近邻搜索
	kdtree->setInputCloud(map); // 设置要搜索的点云，建立KDTree
	std::vector<int> indices;     // 存储查询近邻点索引
	std::vector<float> distances; // 存储近邻点对应距离的平方
	int k = 3;
	vector<pcl::PointCloud<PointT>::Ptr> cloud_vec;
	for (int i = 0; i < MAX_PARTICLES; ++i)
	{
		double theta_tx = particles1[i].x;
		double theta_ty = particles1[i].y;
		double theta_o = particles2[i].o;
		//double theta_o = query.groundtruth.o ;
		double r[2] = { cos(theta_o), -sin(theta_o) };

		//更新当前点云数据
		PointT thisPoint;
		int cloudSize = query.cloud->points.size();
		pcl::PointCloud<PointT>::Ptr cloud_after(new pcl::PointCloud<PointT>);
		double dis = 0;
		for (int j = 0; j < cloudSize; ++j)
		{
			thisPoint.x = (query.cloud->points[j].x) * r[0] - (query.cloud->points[j].y) * r[1] + theta_tx;
			thisPoint.y = (query.cloud->points[j].x) * r[1] + (query.cloud->points[j].y) * r[0] + theta_ty;
			thisPoint.z = 0;
			//thisPoint.intensity = query.cloud->points[j].intensity;
			cloud_after->points.push_back(thisPoint);
			kdtree->nearestKSearch(thisPoint, k, indices, distances);
			// 3.2计算粒子权重
			for (int n = 0; n < k; n++) {
				double x = map->points[indices[n]].x;
				double y = map->points[indices[n]].y;
				double B = (x - thisPoint.x)*(x - thisPoint.x) / (2 * std_x*std_x) + (y - thisPoint.y)*(y - thisPoint.y) / (2 * std_y*std_y);
				dis = dis + exp(-B);
			}
		}
		//visualize_pcd2(map, cloud_after);
		dis_vec.push_back(dis / cloudSize);
		cloud_vec.push_back(cloud_after);
	}
	int maxindex = max_element(dis_vec.begin(), dis_vec.end()) - dis_vec.begin();
	visualize_pcd2(map, cloud_vec[maxindex]);

	/*int minindex = min_element(dis_vec.begin(), dis_vec.end()) - dis_vec.begin();
	visualize_pcd2(map, cloud_vec[minindex]);*/
	double weights_sum1 = 0, weights_sum2 = 0;
	boost::math::normal_distribution<> norm(1, 0.2);
	double weights_dis[MAX_PARTICLES];
	for (int i = 0; i < MAX_PARTICLES; i++)
	{
		weights_dis[i] = pdf(norm, dis_vec[i] / dis_vec[maxindex]);
		particles1[i].w *= weights_dis[i];
		particles2[i].w *= weights_dis[i];
		weights1[i] = particles1[i].w;
		weights2[i] = particles2[i].w;
		weights_sum1 += weights1[i];
		weights_sum2 += weights2[i];
	}
	// normalize weights to bring them in (0, 1]
	for (int i = 0; i < MAX_PARTICLES; i++)
	{
		particles1[i].w /= weights_sum1;
		weights1[i] = particles1[i].w;
		particles2[i].w /= weights_sum2;
		weights2[i] = particles2[i].w;
	}
}