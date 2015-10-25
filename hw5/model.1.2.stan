data{
  vector[2] mu;
  matrix[2,2] cov; 
}
parameters{
  vector[2] params;
}
model {
  params ~ multi_student_t(4,mu,cov);
}