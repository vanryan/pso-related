void individual::fips_update()
{
    int inds = neigh->get_size();
    double coef[inds];
    double xsi = 0.729;
    double phi = 4.01;
    int neighs = 0;

    for(int k = 0; k < inds; k++)
    {
        if(neigh->get_weight(k, index) != 0)
        {
            neighs++;
        }
    }

    for(int d = 0; d < dimensions; d++)
    {
        double sum_coef = 0;
        double pos_ave = 0;

        for(int k = 0; k < inds; k++)
        {
            if(neigh->get_weight(k, index) != 0)
            {
                coef[k] = uniform() * neigh->get_weight(k, index);
                sum_coef += coef[k];
                pos_ave += coef[k] * pop->get_individual(k)->best_position[d];
            }
        }

        pos_ave /= sum_coef;
        double var_phi = sum_coef * phi / neighs;

        double velocity = position[d] - previous_position[d];
        double social_central_tendency = pos_ave - position[d];
        previous_position[d] = position[d];

        position[d] += xsi * (velocity + var_phi * social_central_tendency);
    }
}
